# =============================================================================
# 8b_onset_cessation_dryspells.R
# Onset, cessation, growing season length and dry spells
# Onset / cessation after Dunning, Black & Allan (2016, JGR-Atmospheres,
# doi:10.1002/2016JD025428), annual (monomodal) regime, i.e. the
# cumulative daily rainfall anomaly method of Liebmann et al. (2012).
# Observations vs CHIRPS v2, CHIRPS v3 (RNL), CHIRPS v3 (SAT), 1998-2020
# Outputs: figures/8b_onset_cessation_dryspells.png
#          outputs/8b_season_metrics_station_year.csv
#          outputs/8b_season_climatology.csv
# =============================================================================
#
# Method (applied separately to every station x source series):
#   1. Climatology. Qi = mean rainfall of calendar day i (365 days, 29 Feb
#      dropped); Qbar = mean(Qi), the climatological mean daily rainfall.
#      C(d) = cumsum(Qi - Qbar). ds = DOY of min C (start of the climatological
#      water season), de = DOY of max C (its end).
#   2. Each year. A(D) = cumsum(Rj - Qbar) for j from ds - 50 to de + 50.
#      Onset     = day after the minimum of A(D).
#      Cessation = day of the maximum of A(D), searched after that minimum.
#   Growing season length = cessation - onset (days).
#   Dry spells: runs of days with P < rain_thr between onset and cessation,
#   exactly as in script 8.
# Qbar is computed per series, so each source is judged against its own
# climatology (the method is insensitive to systematic rainfall bias).
# =============================================================================

root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(ggstatsplot)   # robust within-subject stats also need WRS2
  library(patchwork)
})

# ----------------------------- Settings --------------------------------------
years        <- 1998:2020
rain_thr     <- 1      # mm: a day is "wet" if P >= 1 mm (dry spells only)
min_ds       <- 2      # minimum dry spell length (days) reported
win_buffer   <- 50     # days added before ds and after de (Dunning et al., annual regime)
max_na_days  <- 0      # max missing days tolerated in a yearly window; gaps are
# filled with Qbar (zero anomaly, neutral for A(D))
drop_edge    <- TRUE   # discard a year if the min/max of A(D) falls on the
# window boundary (extremum not really resolved)
harm_ratio_max <- 1    # A2/A1 above this = biannual-looking climatology (warning)

f_stations <- "data/csv/p_stations.csv"
f_obs      <- "data/csv/p_stations_daily_data.csv"
chirps_dirs <- c("CHIRPS v2"        = "data/chirps/chirps_v2",
                 "CHIRPS v3 (RNL)"  = "data/chirps/chirps_v3_rnl",
                 "CHIRPS v3 (SAT)"  = "data/chirps/chirps_v3_sat")

src_levels <- c("Observations", names(chirps_dirs))
src_cols   <- c("Observations"    = "grey30",
                "CHIRPS v2"       = "darkred",
                "CHIRPS v3 (RNL)" = "darkblue",
                "CHIRPS v3 (SAT)" = "darkgreen")

all_dates <- seq(as.Date(paste0(min(years), "-01-01")),
                 as.Date(paste0(max(years), "-12-31")), by = "day")

# ----------------------------- Stations --------------------------------------
stn <- read.csv(f_stations, check.names = FALSE, stringsAsFactors = FALSE)
stn$Station <- trimws(stn$Station)

# ----------------------------- Observations ----------------------------------
obs <- read.csv(f_obs, check.names = FALSE, stringsAsFactors = FALSE)
names(obs)[1] <- "Date"
names(obs)[-1] <- trimws(names(obs)[-1])
obs$Date <- as.Date(obs$Date, format = "%m/%d/%Y")
obs <- obs[!is.na(obs$Date) & obs$Date >= min(all_dates) & obs$Date <= max(all_dates), ]

# Keep stations present in both files, same order everywhere
stn_names <- intersect(stn$Station, names(obs)[-1])
if (length(stn_names) < nrow(stn))
  warning("Stations without matching data column: ",
          paste(setdiff(stn$Station, stn_names), collapse = ", "))
stn <- stn[match(stn_names, stn$Station), ]

# Matrix [all_dates x stations]; missing dates become NA
obs_mat <- as.matrix(obs[match(all_dates, obs$Date), stn_names, drop = FALSE])
storage.mode(obs_mat) <- "numeric"
obs_mat[obs_mat < 0] <- NA
rownames(obs_mat) <- as.character(all_dates)

# ----------------------------- CHIRPS extraction -----------------------------
pts <- vect(data.frame(x = stn$Longitude, y = stn$Latitude),
            geom = c("x", "y"), crs = "EPSG:4326")

extract_chirps <- function(dir) {
  files <- list.files(dir, pattern = "[0-9]{4}\\.tif$", full.names = TRUE)
  fy    <- as.integer(sub(".*([0-9]{4})\\.tif$", "\\1", files))
  files <- files[fy %in% years][order(fy[fy %in% years])]
  fy    <- sort(fy[fy %in% years])
  if (!setequal(fy, years)) stop("Missing years in ", dir)
  
  out <- lapply(seq_along(files), function(i) {
    r <- rast(files[i])
    nd <- as.integer(as.Date(paste0(fy[i], "-12-31")) -
                       as.Date(paste0(fy[i], "-01-01"))) + 1L
    if (nlyr(r) != nd)
      stop(basename(files[i]), ": ", nlyr(r), " bands, expected ", nd)
    p <- if (!terra::same.crs(r, pts)) terra::project(pts, terra::crs(r)) else pts
    ex <- as.matrix(terra::extract(r, p, ID = FALSE))   # stations x bands
    t(ex)                                                # bands x stations
  })
  m <- do.call(rbind, out)
  colnames(m) <- stn_names
  rownames(m) <- as.character(all_dates)
  m[m < 0] <- NA                                         # nodata (e.g. -9999)
  m
}

mats <- c(list("Observations" = obs_mat), lapply(chirps_dirs, extract_chirps))

# ----------------------------- Helpers ---------------------------------------
is_leap <- function(y) (y %% 4 == 0 & y %% 100 != 0) | y %% 400 == 0

# Day of year with leap-year correction (so 1 Mar = 60 every year)
doy_norm <- function(d) {
  y <- as.integer(format(d, "%Y")); j <- as.integer(format(d, "%j"))
  j - (is_leap(y) & j > 60)
}

# Inverse of doy_norm on the 365-day calendar (29 Feb never returned)
doy_to_date <- function(y, k) {
  as.Date(sprintf("%d-01-01", y)) + (k - 1) + (is_leap(y) & k >= 60)
}

# ----------------------------- Dunning et al. (2016) --------------------------
# Step 1: climatological water season of one daily series
season_clim <- function(p, dates) {
  keep <- format(dates, "%m-%d") != "02-29"
  Qi   <- tapply(p[keep], doy_norm(dates[keep]), mean, na.rm = TRUE)
  out  <- list(Qbar = NA_real_, ds = NA_integer_, de = NA_integer_,
               harm_ratio = NA_real_)
  if (length(Qi) != 365 || anyNA(Qi)) return(out)
  
  Qbar <- mean(Qi)                      # climatological mean daily rainfall
  C    <- cumsum(Qi - Qbar)             # climatological cumulative anomaly
  f    <- fft(Qi)                       # harmonics: f[2] annual, f[3] semi-annual
  list(Qbar       = Qbar,
       ds         = as.integer(which.min(C)),
       de         = as.integer(which.max(C)),
       harm_ratio = Mod(f[3]) / Mod(f[2]))
}

# Step 2: onset, cessation, length and dry spells for one season (year y)
season_metrics <- function(p, dates, clim, y) {
  res <- c(onset = NA, cess = NA, len = NA, mean_ds = NA, max_ds = NA)
  if (is.na(clim$Qbar)) return(res)
  
  # Window [ds - 50, de + 50]; de may fall in the next year if the
  # climatological season spans 31 December
  y_end <- if (clim$de > clim$ds) y else y + 1L
  i_s <- match(doy_to_date(y,     clim$ds), dates) - win_buffer
  i_e <- match(doy_to_date(y_end, clim$de), dates) + win_buffer
  if (is.na(i_s) || is.na(i_e) || i_s < 1 || i_e > length(p)) return(res)
  
  w  <- i_s:i_e
  pw <- p[w]
  na <- is.na(pw)
  if (sum(na) > max_na_days) return(res)
  pw[na] <- clim$Qbar
  
  A    <- cumsum(pw - clim$Qbar)        # yearly cumulative anomaly A(D)
  n    <- length(A)
  imin <- which.min(A)
  if (imin >= n - 1) return(res)
  imax <- imin + which.max(A[(imin + 1):n])
  if (drop_edge && (imin == 1 || imax == n)) return(res)
  
  onset <- w[imin] + 1L                 # day AFTER the minimum
  cess  <- w[imax]                      # day OF the maximum
  if (cess <= onset) return(res)
  
  res["onset"] <- doy_norm(dates[onset])
  res["cess"]  <- doy_norm(dates[cess])
  res["len"]   <- cess - onset
  
  # Dry spells >= 2 days between onset and cessation (NA days are not dry)
  dry <- !(p[onset:cess] >= rain_thr)
  dry[is.na(dry)] <- FALSE
  r  <- rle(dry)
  ds <- r$lengths[r$values & r$lengths >= min_ds]
  if (length(ds) > 0) {
    res["mean_ds"] <- mean(ds)
    res["max_ds"]  <- max(ds)
  }
  res
}

# ----------------------------- Metrics ---------------------------------------
clim_tab <- list()
metrics <- do.call(rbind, lapply(names(mats), function(s) {
  m <- mats[[s]]
  do.call(rbind, lapply(stn_names, function(st) {
    cl <- season_clim(m[, st], all_dates)
    clim_tab[[length(clim_tab) + 1]] <<- data.frame(
      Source = s, Station = st, Qbar = cl$Qbar, ds = cl$ds, de = cl$de,
      harm_ratio = cl$harm_ratio)
    do.call(rbind, lapply(years, function(y) {
      data.frame(Source = s, Station = st, Year = y,
                 t(season_metrics(m[, st], all_dates, cl, y)),
                 row.names = NULL)
    }))
  }))
}))
metrics$Source <- factor(metrics$Source, levels = src_levels)
clim_tab <- do.call(rbind, clim_tab)

bimodal <- clim_tab[!is.na(clim_tab$harm_ratio) & clim_tab$harm_ratio > harm_ratio_max, ]
if (nrow(bimodal) > 0)
  warning("Harmonic ratio A2/A1 > ", harm_ratio_max,
          " (not clearly monomodal) for: ",
          paste(bimodal$Source, bimodal$Station, sep = " / ", collapse = "; "))

dir.create("outputs", showWarnings = FALSE)
write.csv(metrics,  "outputs/8b_season_metrics_station_year.csv", row.names = FALSE)
write.csv(clim_tab, "outputs/8b_season_climatology.csv",         row.names = FALSE)

# Per-year values: mean across stations (one value per source and year)
yearly <- aggregate(cbind(onset, cess, len, mean_ds, max_ds) ~ Source + Year,
                    data = metrics, FUN = mean, na.rm = TRUE, na.action = na.pass)
for (v in c("onset", "cess", "len", "mean_ds", "max_ds"))
  yearly[[v]][is.nan(yearly[[v]])] <- NA   # all-NA years -> NA instead of NaN

# ----------------------------- Plot ------------------------------------------
txt_theme <- theme(
  plot.title = element_text(size = 18, colour = "black", hjust = 0, face = "bold"),
  axis.title = element_text(colour = "black", size = 12),
  axis.text  = element_text(colour = "black", size = 12)
)

make_panel <- function(var, title, ylab, show_x = TRUE) {
  d <- yearly[, c("Source", "Year", var)]
  names(d)[3] <- "value"
  # keep only years available for all four sources (paired design)
  ok <- tapply(!is.na(d$value), d$Year, sum) == length(src_levels)
  d  <- d[d$Year %in% as.integer(names(ok)[ok]), ]
  d  <- d[order(d$Source, d$Year), ]
  d$Source <- droplevels(d$Source)
  
  p <- ggwithinstats(
    data = d, x = Source, y = value,
    type = "robust",
    pairwise.display = "none",
    bf.message = FALSE,
    title = title, xlab = NULL, ylab = ylab
  ) +
    scale_color_manual(values = src_cols) +
    scale_fill_manual(values = src_cols) +
    txt_theme +
    theme(axis.title.x = element_blank())          # no "Source" title anywhere
  
  if (!show_x)                                     # x text only on bottom-most panels
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  p
}

# Layout is 2 x 3 with an empty slot (row 2, col 3): the bottom-most panels are
# d (col 1), e (col 2) and c (col 3); only a and b hide their x text.
p_a <- make_panel("onset",   "a) Onset date",               "DOY", show_x = FALSE)
p_b <- make_panel("cess",    "b) Cessation date",           "DOY", show_x = FALSE)
p_c <- make_panel("len",     "c) Growing season length",    "Days")
p_d <- make_panel("mean_ds", "d) Dry spell duration (avg)", "Days")
p_e <- make_panel("max_ds",  "e) Longest dry spell (avg)",  "Days")

fig <- wrap_plots(p_a, p_b, p_c, p_d, p_e, nrow = 2, ncol = 3)

dir.create("figures", showWarnings = FALSE)
ggsave("figures/8_onset_cessation_dryspells.png", fig,
       width = 25, height = 14, units = "in",
       dpi = 500, bg = "white", scale = 0.8)
