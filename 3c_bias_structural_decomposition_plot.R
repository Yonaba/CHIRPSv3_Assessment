root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# ------------------------------------------------------------------
# Decomposition of hit biases (over / under) of CHIRPS v2, v3 (RNL)
# and v3 (SAT) into systematic and random components, by month and
# by year. Based on a log-log regression: log(CHIRPS) ~ log(OBS)
# ------------------------------------------------------------------

library(terra)
library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

# ---- Settings ----------------------------------------------------
period   <- 1998:2020
MIN_TH   <- 0.1          # rain / no-rain threshold (mm)
MONTHS   <- 6:10         # months analysed in the monthly figure (Jun-Oct)
MIN_N    <- 3            # minimum number of hits needed to fit a regression
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
stations   <- stations[match(st_names, stations$Station), ]
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

# ---- 3. Decomposition of hit biases ------------------------------
# Dates / months / years aligned with the vectorised [days x stations] matrices
dates_vec <- rep(all_dates, times = N_STATIONS)
mon_vec   <- as.numeric(format(dates_vec, "%m"))
year_vec  <- as.numeric(format(dates_vec, "%Y"))

# Systematic / random shares (%) of the squared log error, for one subset
# of hits (over- or under-estimations) of one group (month or year)
decomp <- function(d) {
  if (nrow(d) < MIN_N || length(unique(d$OBS)) < 2) return(c(sys = NA, rnd = NA))
  fit  <- lm(log(SIM) ~ log(OBS), data = d)
  rlin <- exp(as.numeric(predict(fit, data.frame(OBS = d$OBS))))
  sys  <- mean(log(rlin / d$OBS)^2, na.rm = TRUE)
  rnd  <- mean(log(d$SIM / rlin)^2, na.rm = TRUE)
  tot  <- sys + rnd
  c(sys = sys / tot * 100, rnd = rnd / tot * 100)
}

decompose_by <- function(group_vec, groups) {
  res <- list()
  for (p in names(products)) {
    message("Decomposing biases: ", p)
    d <- data.frame(grp = group_vec,
                    OBS = round(as.vector(obs_mat), 1),
                    SIM = round(as.vector(chirps[[p]]), 1))
    d <- d[complete.cases(d), ]
    # hits with different amounts
    d <- d[d$OBS >= MIN_TH & d$SIM >= MIN_TH & d$SIM != d$OBS, ]
    for (g in groups) {
      dg <- d[d$grp == g, ]
      over  <- decomp(dg[dg$SIM > dg$OBS, ])
      under <- decomp(dg[dg$SIM < dg$OBS, ])
      res[[length(res) + 1]] <- data.frame(
        product = p, group = g,
        bias_type = c("Over Hit Bias (OHB)", "Under Hit Bias (UHB)"),
        sys = c(over[["sys"]],  -under[["sys"]]),   # under-estimation shown as negative
        rnd = c(over[["rnd"]],  -under[["rnd"]]))
    }
  }
  bind_rows(res)
}

bias_mon <- decompose_by(mon_vec,  MONTHS)
bias_ann <- decompose_by(year_vec, period)

write_csv(bias_mon, "tables/2_hit_bias_decomposition_monthly.csv")
write_csv(bias_ann, "tables/2_hit_bias_decomposition_annual.csv")

# ---- 4. Plots -----------------------------------------------------
fill_levels <- c("Over Hit Bias (OHB) - Random",
                 "Under Hit Bias (UHB) - Random",
                 "Over Hit Bias (OHB) - Systematic",
                 "Under Hit Bias (UHB) - Systematic")
fill_values <- setNames(c("blue3", "red3", "#95D0FC", "orange"), fill_levels)

make_plots <- function(bias_df, x_breaks, x_labels, rotate = FALSE) {
  plist <- list()
  for (cp in seq_along(products)) {
    p <- names(products)[cp]
    df <- bias_df[bias_df$product == p, ]
    long <- bind_rows(
      transform(df[, c("group", "bias_type")], value = df$sys, err = "Systematic"),
      transform(df[, c("group", "bias_type")], value = df$rnd, err = "Random"))
    long$fill <- factor(paste(long$bias_type, long$err, sep = " - "),
                        levels = fill_levels)
    # stacking order: random first, then systematic
    long <- long[order(long$fill), ]
    
    plist[[cp]] <- ggplot(long, aes(x = group, y = value, fill = fill)) +
      geom_area() +
      ylab(ifelse(cp == 1, "Percentage (%)", "")) +
      xlab("") +
      scale_x_continuous(breaks = x_breaks, labels = x_labels) +
      scale_fill_manual(name = "Hit bias decomposition (%)",
                        values = fill_values, drop = FALSE) +
      labs(title = paste0("(", letters[cp], ") ", p)) +
      theme_bw() +
      theme(plot.title = element_text(color = "black", size = 18, face = "bold"),
            axis.title.y = element_text(size = 12, margin = margin(r = 10)),
            axis.title.x = element_text(size = 12, margin = margin(t = 10)),
            axis.text = element_text(color = "black", size = 12),
            axis.text.x = if (rotate) element_text(angle = 90, vjust = 0.5) else element_text(),
            legend.position = "bottom",
            legend.title = element_text(size = 12, color = "black"),
            legend.text = element_text(size = 12, color = "black")) +
      guides(fill = guide_legend(nrow = 2))
  }
  ggarrange(plotlist = plist, nrow = 1, ncol = 3,
            common.legend = TRUE, legend = "bottom")
}

grob_mon <- make_plots(bias_mon, MONTHS, month.abb[MONTHS])
ggsave("figures/3_hit_bias_decomposition_months.png", plot = grob_mon,
       width = 20, height = 7.5, units = "in", dpi = 350, bg = "white", scale = 0.8)

grob_ann <- make_plots(bias_ann, period, period, rotate = TRUE)
ggsave("figures/3_hit_bias_decomposition_annual.png", plot = grob_ann,
       width = 26, height = 7.5, units = "in", dpi = 350, bg = "white", scale = 0.8)
