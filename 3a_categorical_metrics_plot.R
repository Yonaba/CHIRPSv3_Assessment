root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# ------------------------------------------------------------------
# Categorical scores (POD, FAR, TS) of CHIRPS v2, v3 (RNL) and v3 (SAT)
# against observed daily rainfall, for various rainfall thresholds
# ------------------------------------------------------------------

library(terra)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggpubr)

# ---- Settings ----------------------------------------------------
period    <- 1998:2020
products  <- c("CHIRPS v2"      = "chirps_v2",
               "CHIRPS v3 (RNL)" = "chirps_v3_rnl",
               "CHIRPS v3 (SAT)" = "chirps_v3_sat")
cols      <- c("CHIRPS v2"       = "darkred",
               "CHIRPS v3 (RNL)" = "darkblue",
               "CHIRPS v3 (SAT)" = "darkgreen")
# Rain/no-rain threshold used to clean CHIRPS values
min_rain  <- 0.1
# Thresholds (mm/day). An event is "rain >= threshold", except for the last
# class (">100"), where an event is "rain > 100".
thresholds <- c(0.1, 5, 10, 20, 30, 50, 75, 100, 100)
thr_labels <- c("0.1", "5", "10", "20", "30", "50", "75", "100", ">100")
thr_strict <- c(rep(FALSE, 8), TRUE)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load stations and observations ---------------------------
stations <- read_csv("data/csv/p_stations.csv", show_col_types = FALSE) %>%
  as.data.frame()

obs_wide <- read_csv("data/csv/p_stations_daily_data.csv",
                     show_col_types = FALSE, name_repair = "minimal") %>%
  as.data.frame()
names(obs_wide)[1] <- "Date"
obs_wide$Date <- as.Date(obs_wide$Date, format = "%m/%d/%Y")
obs_wide <- obs_wide %>%
  filter(Date >= as.Date("1998-01-01"), Date <= as.Date("2020-12-31"))

# Keep only stations present in both files, in the same order
st_names <- intersect(stations$Station, names(obs_wide)[-1])
if (length(st_names) < nrow(stations)) {
  warning("Stations without coordinates or observations were dropped: ",
          paste(setdiff(stations$Station, st_names), collapse = ", "))
}
stations <- stations[match(st_names, stations$Station), ]

all_dates <- seq(as.Date("1998-01-01"), as.Date("2020-12-31"), by = "day")
obs_wide  <- obs_wide[match(all_dates, obs_wide$Date), c("Date", st_names)]

# Observations: matrix [days x stations]
obs_mat <- as.matrix(obs_wide[, st_names])
storage.mode(obs_mat) <- "numeric"

# ---- 2. Extract CHIRPS time series at stations -------------------
pts <- vect(stations, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")

extract_chirps <- function(folder) {
  files <- list.files(file.path("data/chirps", folder),
                      pattern = "[0-9]{4}\\.tif$", full.names = TRUE)
  years <- as.integer(sub(".*([0-9]{4})\\.tif$", "\\1", files))
  files <- files[years %in% period][order(years[years %in% period])]
  years <- sort(years[years %in% period])
  
  out <- lapply(seq_along(files), function(i) {
    r <- rast(files[i])
    dts <- seq(as.Date(paste0(years[i], "-01-01")), by = "day",
               length.out = nlyr(r))
    if (nlyr(r) != ifelse(leap_year(years[i]), 366, 365)) {
      warning(basename(files[i]), ": unexpected number of bands (", nlyr(r), ")")
    }
    vals <- as.matrix(terra::extract(r, pts, ID = FALSE))   # stations x days
    vals <- t(vals)                                  # days x stations
    rownames(vals) <- as.character(dts)
    vals
  })
  mat <- do.call(rbind, out)
  colnames(mat) <- stations$Station
  mat <- mat[match(as.character(all_dates), rownames(mat)), , drop = FALSE]
  # Fill value / negative values -> NA, then values < 0.1 -> 0
  mat[!is.na(mat) & mat < 0] <- NA
  mat[!is.na(mat) & mat < min_rain] <- 0
  mat
}

leap_year <- function(y) (y %% 4 == 0 & y %% 100 != 0) | y %% 400 == 0

chirps <- lapply(products, extract_chirps)
names(chirps) <- names(products)

# ---- 3. Categorical scores ---------------------------------------
compute_scores <- function(sim_mat, obs_mat, thr, strict) {
  ok  <- !is.na(sim_mat) & !is.na(obs_mat)
  sim <- sim_mat[ok]
  obs <- obs_mat[ok]
  if (strict) { o <- obs > thr; s <- sim > thr } else { o <- obs >= thr; s <- sim >= thr }
  # Contingency table counts, stored as doubles: verification::table.stats()
  # uses integer counts and overflows with large samples (stations x days).
  a <- as.numeric(sum(o & s))    # hits
  b <- as.numeric(sum(!o & s))   # false alarms
  c <- as.numeric(sum(o & !s))   # misses
  d <- as.numeric(sum(!o & !s))  # correct negatives
  # Same fudge as verification::table.stats(): if any cell is empty,
  # a small number is added to every cell
  fudge <- 0.01
  if (any(c(a, b, c, d) == 0)) { a <- a + fudge; b <- b + fudge
  c <- c + fudge; d <- d + fudge }
  data.frame(POD = a / (a + c), FAR = b / (a + b), TS = a / (a + b + c))
}

scores <- bind_rows(lapply(names(chirps), function(p) {
  bind_rows(lapply(seq_along(thresholds), function(k) {
    cbind(Product = p,
          Threshold = thr_labels[k],
          compute_scores(chirps[[p]], obs_mat, thresholds[k], thr_strict[k]))
  }))
}))

scores$Product   <- factor(scores$Product, levels = names(products))
scores$Threshold <- factor(scores$Threshold, levels = thr_labels)

write_csv(scores, "tables/2_scores_chirps_vs_obs.csv")

# ---- 4. Plot ------------------------------------------------------
scores_long <- scores %>%
  pivot_longer(c(POD, FAR, TS), names_to = "Metric", values_to = "Score")

metric_titles <- c(POD = "a) Probability of detection (POD)",
                   FAR = "b) False alarm ratio (FAR)",
                   TS  = "c) Threat score (TS)")

make_panel <- function(m) {
  ggplot(filter(scores_long, Metric == m),
         aes(x = Threshold, y = Score, colour = Product, group = Product)) +
    geom_line(linetype = "dashed", linewidth = 0.8) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = metric_titles[[m]],
         x = "\nRainfall threshold (mm/day)", y = ifelse(m=="POD", "Score\n", "")) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 18, hjust = 0, colour = "black"),
          axis.title = element_text(size = 12, colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          legend.title = element_text(size = 12, colour = "black"),
          legend.text = element_text(size = 12, colour = "black"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
}

panels <- lapply(c("POD", "FAR", "TS"), make_panel)

fig <- ggarrange(plotlist = panels, ncol = 3, nrow = 1,
                 common.legend = TRUE, legend = "bottom")

ggsave("figures/3_chirps_categorical_scores.png", fig,
       width = 15, height = 5.5, units = "in", dpi = 350, bg = "white", scale = 0.95)
