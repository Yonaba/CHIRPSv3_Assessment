root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# ------------------------------------------------------------------
# Hit / miss / false-alarm biases of CHIRPS v2, v3 (RNL) and v3 (SAT)
# against observed daily rainfall (event and annual average biases)
# ------------------------------------------------------------------

library(terra)
library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

# ---- Settings ----------------------------------------------------
period   <- 1998:2020
N_YEARS  <- length(period)
MIN_TH   <- 0.1     # rain / no-rain threshold (mm)
products <- c("CHIRPS v2"       = "chirps_v2",
              "CHIRPS v3 (RNL)" = "chirps_v3_rnl",
              "CHIRPS v3 (SAT)" = "chirps_v3_sat")

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

leap_year <- function(y) (y %% 4 == 0 & y %% 100 != 0) | y %% 400 == 0

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

st_names <- intersect(stations$Station, names(obs_wide)[-1])
if (length(st_names) < nrow(stations)) {
  warning("Stations without coordinates or observations were dropped: ",
          paste(setdiff(stations$Station, st_names), collapse = ", "))
}
stations <- stations[match(st_names, stations$Station), ]
N_STATIONS <- nrow(stations)

all_dates <- seq(as.Date("1998-01-01"), as.Date("2020-12-31"), by = "day")
obs_wide  <- obs_wide[match(all_dates, obs_wide$Date), c("Date", st_names)]
obs_mat   <- as.matrix(obs_wide[, st_names])
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
    vals <- t(as.matrix(terra::extract(r, pts, ID = FALSE)))  # days x stations
    rownames(vals) <- as.character(dts)
    vals
  })
  mat <- do.call(rbind, out)
  colnames(mat) <- stations$Station
  mat <- mat[match(as.character(all_dates), rownames(mat)), , drop = FALSE]
  mat[!is.na(mat) & mat < 0] <- NA          # fill values -> NA
  mat[!is.na(mat) & mat < MIN_TH] <- 0      # values < 0.1 -> 0
  mat
}

chirps <- lapply(products, extract_chirps)
names(chirps) <- names(products)

# ---- 3. Biases ----------------------------------------------------
# Values are rounded to 0.1 mm. All stations and days are pooled.
# "Average" = mean bias per event; "volume" = average annual bias per
# station (total bias / (number of stations x number of years)).
ffdiv <- N_STATIONS * N_YEARS

bias_list <- list()
for (p in names(products)) {
  message("Processing biases: ", p)
  
  d <- data.frame(OBS = round(as.vector(obs_mat), 1),
                  SIM = round(as.vector(chirps[[p]]), 1))
  d <- d[complete.cases(d), ]
  n_all <- nrow(d)
  
  # Hit bias (both detect rain, different amounts)
  h  <- d[d$OBS >= MIN_TH & d$SIM >= MIN_TH & d$SIM != d$OBS, ]
  hbias <- sum(abs(h$SIM - h$OBS))
  
  oh <- h[h$SIM > h$OBS, ]                     # over-estimated hits
  uh <- h[h$SIM < h$OBS, ]                     # under-estimated hits
  ohbias <- sum(oh$SIM - oh$OBS)
  uhbias <- sum(uh$OBS - uh$SIM)
  
  # Miss bias (rain observed, none in CHIRPS)
  m <- d[d$OBS >= MIN_TH & d$SIM < MIN_TH, ]
  mbias <- sum(abs(m$SIM - m$OBS))
  
  # False alarm bias (no observed rain, rain in CHIRPS)
  f <- d[d$OBS < MIN_TH & d$SIM >= MIN_TH, ]
  fabias <- sum(f$SIM - f$OBS)
  
  tbias <- hbias + mbias + fabias
  
  bias_list[[p]] <- data.frame(
    product = p,
    bias_type = c("Hit Bias", "Over Hit Bias", "Under Hit Bias",
                  "Miss Bias", "False Bias", "Total Bias"),
    total   = c(hbias, ohbias, uhbias, mbias, fabias, tbias),
    n       = c(nrow(h), nrow(oh), nrow(uh), nrow(m), nrow(f), n_all)
  )
}

bias.df <- bind_rows(bias_list) %>%
  mutate(average = ifelse(n > 0, total / n, 0),
         volume  = total / ffdiv)

btypes <- c("Hit Bias", "Over Hit Bias", "Under Hit Bias",
            "Miss Bias", "False Bias", "Total Bias")
bias.df$product   <- ordered(bias.df$product, levels = names(products))
bias.df$bias_type <- ordered(bias.df$bias_type, levels = btypes)
bias.df$average   <- round(bias.df$average, 1)
bias.df$volume    <- round(bias.df$volume, 1)

write_csv(bias.df, "figures/3_biases_chirps_vs_obs.csv")

# ---- 4. Plot ------------------------------------------------------
acro <- c("Hit Bias"        = "(HB)",
          "Over Hit Bias"   = "(OHB)",
          "Under Hit Bias"  = "(-UHB)",
          "Miss Bias"       = "(-MB)",
          "False Bias"      = "(FB)",
          "Total Bias"      = "(TB)")

# Scaling of the secondary axis (annual bias) onto the primary axis
coeff <- max(10, ceiling(max(bias.df$volume) / max(bias.df$average) / 10) * 10)
ymax  <- ceiling(max(bias.df$average, bias.df$volume / coeff) / 5) * 5

ctheme <- theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
                axis.title.x = element_text(size = 14, margin = margin(t = 10)),
                axis.title.y = element_text(size = 14, margin = margin(r = 10), colour = "blue"),
                axis.title.y.right = element_text(size = 14, margin = margin(l = 10), colour = "red"),
                axis.line.y = element_line(colour = "blue"),
                axis.line.y.right = element_line(colour = "red"),
                axis.text = element_text(color = "black", size = 12),
                axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 10, colour = "black"),
                axis.text.y = element_text(colour = "blue"),
                axis.text.y.right = element_text(colour = "red"),
                axis.ticks.y.left = element_line(colour = "blue"),
                axis.ticks.y.right = element_line(colour = "red"))

plist <- list()
for (cp in seq_along(btypes)) {
  btype <- btypes[cp]
  ddf <- bias.df[bias.df$bias_type == btype, ]
  plist[[cp]] <- ggplot(data = ddf, aes(x = product)) +
    geom_bar(aes(y = average), fill = "blue", width = 0.4, stat = "identity",
             just = 1, colour = "black", linewidth = 0.5) +
    geom_bar(aes(y = volume / coeff), fill = "red", width = 0.4, stat = "identity",
             just = 0, colour = "black", linewidth = 0.5) +
    xlab("") +
    scale_y_continuous(
      name = ifelse(cp %in% c(1, 4), "Avg. event bias [mm]", ""),
      limits = c(0, ymax),
      sec.axis = sec_axis(~ . * coeff,
                          name = ifelse(cp %in% c(3, 6), "Avg. annual bias [mm]", ""))
    ) +
    labs(title = bquote(bold("(" * .(letters[cp]) * ")" ~ .(btype) ~
                               bolditalic(.(acro[[btype]]))))) +
    theme_bw() + ctheme
}

grob <- ggarrange(plotlist = plist, nrow = 2, ncol = 3)

ggsave("figures/3_chirps_hit_miss_false_biases.png", plot = grob,
       width = 22, height = 12, dpi = 350, scale = 0.75, bg = "white")

