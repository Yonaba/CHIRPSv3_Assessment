# =============================================================================
# 8_onset_cessation_dryspells.R
# Onset, cessation, growing season length and dry spells (Sivakumar, 1988)
# Observations vs CHIRPS v2, CHIRPS v3 (RNL), CHIRPS v3 (SAT), 1998-2020
# Output: figures/8_onset_cessation_dryspells.png
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
rain_thr     <- 1      # mm: a day is "wet" if P >= 1 mm
onset_acc    <- 20     # mm over 3 consecutive days
onset_win    <- 30     # days checked after the onset day for dry spells
onset_maxds  <- 7      # max tolerated dry spell (days) in that window
cess_win     <- 20     # days without rain after the cessation date
min_ds       <- 2      # minimum dry spell length (days) reported

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
    t(ex)                                        # bands x stations
  })
  m <- do.call(rbind, out)
  colnames(m) <- stn_names
  rownames(m) <- as.character(all_dates)
  m[m < 0] <- NA                                 # nodata (e.g. -9999)
  m
}

mats <- c(list("Observations" = obs_mat), lapply(chirps_dirs, extract_chirps))

# ----------------------------- Metrics ---------------------------------------
# Day of year with leap-year correction (so 1 Mar = 60 every year)
doy_norm <- function(d) {
  y <- as.integer(format(d, "%Y")); j <- as.integer(format(d, "%j"))
  leap <- (y %% 4 == 0 & y %% 100 != 0) | y %% 400 == 0
  j - (leap & j > 60)
}

max_dry_run <- function(wet) {      # longest run of dry days in a logical vector
  r <- rle(!wet)
  if (any(r$values)) max(r$lengths[r$values]) else 0L
}

season_metrics <- function(p, dates) {
  n    <- length(p)
  wet  <- p >= rain_thr
  may1 <- which(format(dates, "%m-%d") == "05-01")
  sep1 <- which(format(dates, "%m-%d") == "09-01")
  res  <- c(onset = NA, cess = NA, len = NA, mean_ds = NA, max_ds = NA)
  
  # Onset: first day after 1 May with >= 20 mm over 3 consecutive days
  # (that day and the next two) and no dry spell > 7 days in the next 30 days
  onset <- NA_integer_
  for (d in may1:(n - onset_win)) {
    p3 <- p[d:(d + 2)]
    if (anyNA(p3) || sum(p3) < onset_acc) next
    w <- wet[(d + 1):(d + onset_win)]
    if (anyNA(w)) next
    if (max_dry_run(w) <= onset_maxds) { onset <- d; break }
  }
  
  # Cessation: first day on/after 1 Sep followed by 20 days without rain
  cess <- NA_integer_
  for (d in sep1:(n - cess_win)) {
    w <- wet[(d + 1):(d + cess_win)]
    if (anyNA(w)) next
    if (!any(w)) { cess <- d; break }
  }
  
  if (is.na(onset) || is.na(cess) || cess <= onset) return(res)
  
  res["onset"] <- doy_norm(dates[onset])
  res["cess"]  <- doy_norm(dates[cess])
  res["len"]   <- cess - onset
  
  # Dry spells >= 2 days between onset and cessation (NA days are not dry)
  dry <- !wet[onset:cess]
  dry[is.na(dry)] <- FALSE
  r  <- rle(dry)
  ds <- r$lengths[r$values & r$lengths >= min_ds]
  if (length(ds) > 0) {
    res["mean_ds"] <- mean(ds)
    res["max_ds"]  <- max(ds)
  }
  res
}

year_idx <- split(seq_along(all_dates), format(all_dates, "%Y"))

metrics <- do.call(rbind, lapply(names(mats), function(s) {
  m <- mats[[s]]
  do.call(rbind, lapply(stn_names, function(st) {
    do.call(rbind, lapply(names(year_idx), function(y) {
      i <- year_idx[[y]]
      data.frame(Source = s, Station = st, Year = as.integer(y),
                 t(season_metrics(m[i, st], all_dates[i])),
                 row.names = NULL)
    }))
  }))
}))
metrics$Source <- factor(metrics$Source, levels = src_levels)

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
p_a <- make_panel("onset",   "a) Onset date",              "DOY", show_x = FALSE)
p_b <- make_panel("cess",    "b) Cessation date",          "DOY", show_x = FALSE)
p_c <- make_panel("len",     "c) Growing season length",   "Days")
p_d <- make_panel("mean_ds", "d) Dry spell duration (avg)", "Days")
p_e <- make_panel("max_ds",  "e) Longest dry spell (avg)",  "Days")

fig <- wrap_plots(p_a, p_b, p_c, p_d, p_e, nrow = 2, ncol = 3)

dir.create("figures", showWarnings = FALSE)
ggsave("figures/8_onset_cessation_dryspells.png", fig,
       width = 24, height = 14, units = "in",
       dpi = 350, bg = "white", scale = 0.8)