root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# Evaluation of CHIRPS v2, CHIRPS v3 (RNL) and CHIRPS v3 (SAT) against station
# observations at 5 timescales (daily, pentad, dekad, month, year), 1998-2020.
# Metrics per station: Pearson r, RMSE, PBIAS, KGE (Kling et al. 2012).
# Output: multi-panel figure (5 rows = timescales x 4 columns = metrics)
#         in "figures/2_*.png"
# =============================================================================

# ---- 0. Packages -------------------------------------------------------------
pkgs <- c("terra", "dplyr", "tidyr", "lubridate", "ggplot2",
          "ggstatsplot", "ggpubr", "readr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 1. Settings -------------------------------------------------------------
start_date <- as.Date("1998-01-01")
end_date   <- as.Date("2020-12-31")
all_dates  <- seq(start_date, end_date, by = "day")

stations_file <- "data/csv/p_stations.csv"
obs_file      <- "data/csv/p_stations_daily_data.csv"
chirps_dir    <- "data/chirps"

products <- c("CHIRPS v2"       = "chirps_v2",
              "CHIRPS v3 (RNL)" = "chirps_v3_rnl",
              "CHIRPS v3 (SAT)" = "chirps_v3_sat")

timescales <- c("Daily", "Pentad", "Dekadal", "Monthly", "Annual")
ts_unit    <- c(Daily = "mm/day", Pentad = "mm/pentad", Dekadal = "mm/dekad",
                Monthly = "mm/month", Annual = "mm/year")

# Minimum number of paired (complete) values required to compute the metrics
# of a station at a given timescale (about one year of data; 72 pentads and
# 36 dekads per calendar year)
min_pairs <- c(Daily = 365, Pentad = 72, Dekadal = 36, Monthly = 24, Annual = 8)

dir.create("figures", showWarnings = FALSE)

# ---- 2. Load station locations ----------------------------------------------
stations <- read_csv(stations_file, show_col_types = FALSE) |>
  mutate(Station = trimws(as.character(Station))) |>
  filter(!is.na(Latitude), !is.na(Longitude))

# ---- 3. Load observed daily precipitation (wide -> long) ---------------------
obs_wide <- read_csv(obs_file, show_col_types = FALSE, name_repair = "minimal")
names(obs_wide) <- trimws(names(obs_wide))
names(obs_wide)[1] <- "Date"

obs_long <- obs_wide |>
  mutate(Date = as.Date(as.character(Date), format = "%m/%d/%Y")) |>
  filter(Date >= start_date, Date <= end_date) |>
  mutate(across(-Date, as.numeric)) |>
  pivot_longer(-Date, names_to = "Station", values_to = "obs")

# Keep stations present in both tables
common_st <- intersect(stations$Station, unique(obs_long$Station))
if (length(common_st) == 0) stop("No station names match between the two csv files.")
if (length(common_st) < nrow(stations))
  message(nrow(stations) - length(common_st),
          " station(s) from p_stations.csv have no data column and are dropped.")
stations <- stations |> filter(Station %in% common_st)

# Complete daily grid (missing days -> NA, so that accumulations can be validated)
obs_long <- tidyr::expand_grid(Date = all_dates, Station = stations$Station) |>
  left_join(obs_long, by = c("Date", "Station"))

# ---- 4. Extract CHIRPS daily series at station locations --------------------
extract_chirps <- function(folder) {
  files <- list.files(file.path(chirps_dir, folder), pattern = "[0-9]{4}\\.tif$",
                      full.names = TRUE)
  yrs <- as.integer(sub(".*([0-9]{4})\\.tif$", "\\1", files))
  keep <- yrs >= year(start_date) & yrs <= year(end_date)
  files <- files[keep]; yrs <- yrs[keep]
  files <- files[order(yrs)];  yrs <- sort(yrs)
  if (length(files) == 0) stop("No tif files found in ", folder)
  
  xy <- as.matrix(stations[, c("Longitude", "Latitude")])
  res <- lapply(seq_along(files), function(i) {
    r <- rast(files[i])
    ndays <- as.integer(as.Date(paste0(yrs[i], "-12-31")) -
                          as.Date(paste0(yrs[i], "-01-01"))) + 1
    if (nlyr(r) != ndays)
      stop("File ", basename(files[i]), " has ", nlyr(r), " bands, expected ", ndays)
    pts <- vect(xy, crs = "EPSG:4326")
    if (!is.lonlat(r)) pts <- project(pts, crs(r))
    vals <- terra::extract(r, pts)                       # stations x days
    vals <- as.matrix(vals[, setdiff(names(vals), "ID"), drop = FALSE])
    vals[vals < 0] <- NA                                 # NoData (-9999)
    dts <- seq(as.Date(paste0(yrs[i], "-01-01")), by = "day", length.out = ndays)
    data.frame(Date = rep(dts, each = nrow(vals)),
               Station = rep(stations$Station, times = ndays),
               sim = as.vector(vals))
  })
  bind_rows(res)
}

sim_list <- lapply(products, function(f) {
  message("Extracting ", f, " ...")
  extract_chirps(f)
})

# ---- 5. Metrics --------------------------------------------------------------
# Kling et al. (2012) KGE; PBIAS > 0 means overestimation by the product
calc_metrics <- function(obs, sim) {
  ok <- is.finite(obs) & is.finite(sim)
  o <- obs[ok]; s <- sim[ok]
  out <- c(r = NA_real_, RMSE = NA_real_, PBIAS = NA_real_, KGE = NA_real_)
  if (length(o) < 3) return(out)
  out["RMSE"]  <- sqrt(mean((s - o)^2))
  out["PBIAS"] <- if (sum(o) != 0) 100 * sum(s - o) / sum(o) else NA_real_
  if (sd(o) > 0 && sd(s) > 0 && mean(o) > 0 && mean(s) > 0) {
    r     <- cor(o, s, method = "pearson")
    beta  <- mean(s) / mean(o)
    gamma <- (sd(s) / mean(s)) / (sd(o) / mean(o))
    out["r"]   <- r
    out["KGE"] <- 1 - sqrt((r - 1)^2 + (beta - 1)^2 + (gamma - 1)^2)
  }
  out
}

# Accumulation key for each timescale (CHIRPS calendar-month based periods)
# - Pentad : 6 per month (1-5, 6-10, 11-15, 16-20, 21-25, 26-end of month)
# - Dekadal: 3 per month (1-10, 11-20, 21-end of month)
add_period <- function(df, scale) {
  df <- switch(scale,
               Daily   = df |> mutate(period = as.character(Date)),
               Pentad  = df |> mutate(period = paste(year(Date), month(Date),
                                                     pmin(ceiling(mday(Date) / 5), 6))),
               Dekadal = df |> mutate(period = paste(year(Date), month(Date),
                                                     pmin(ceiling(mday(Date) / 10), 3))),
               Monthly = df |> mutate(period = paste(year(Date), month(Date))),
               Annual  = df |> mutate(period = as.character(year(Date))))
  df
}

# Accumulate: a period is valid only if all its days are valid (paired obs & sim)
accumulate <- function(df, scale) {
  if (scale == "Daily") return(df |> select(Station, obs, sim))
  df |>
    add_period(scale) |>
    group_by(Station, period) |>
    summarise(obs = if (anyNA(obs)) NA_real_ else sum(obs),
              sim = if (anyNA(sim)) NA_real_ else sum(sim),
              .groups = "drop")
}

metrics_all <- list()
for (p in names(products)) {
  message("Computing metrics for ", p, " ...")
  paired <- obs_long |>
    inner_join(sim_list[[p]], by = c("Date", "Station")) |>
    # same sample for obs & sim: a day is valid only if both are available
    mutate(bad = is.na(obs) | is.na(sim),
           obs = ifelse(bad, NA_real_, obs),
           sim = ifelse(bad, NA_real_, sim)) |>
    select(-bad)
  
  for (sc in timescales) {
    agg <- accumulate(paired, sc)
    m <- agg |>
      group_by(Station) |>
      group_modify(~ {
        n <- sum(is.finite(.x$obs) & is.finite(.x$sim))
        res <- if (n >= min_pairs[[sc]]) calc_metrics(.x$obs, .x$sim) else
          c(r = NA_real_, RMSE = NA_real_, PBIAS = NA_real_, KGE = NA_real_)
        as_tibble(as.list(res))
      }) |>
      ungroup() |>
      mutate(Source = p, Timescale = sc)
    metrics_all[[paste(p, sc)]] <- m
  }
}

metrics <- bind_rows(metrics_all) |>
  pivot_longer(c(r, RMSE, PBIAS, KGE), names_to = "Metric", values_to = "value") |>
  mutate(Source    = factor(Source, levels = names(products)),
         Timescale = factor(Timescale, levels = timescales),
         Metric    = factor(Metric, levels = c("r", "RMSE", "PBIAS", "KGE")))

write_csv(metrics, "data/csv/2_station_metrics_long.csv")

# ---- 6. Figure ---------------------------------------------------------------
metric_names <- c(r = "Pearson r", RMSE = "RMSE", PBIAS = "PBIAS", KGE = "KGE")
ylab_fun <- function(m, sc) {
  switch(as.character(m),
         r     = "r (-)",
         RMSE  = paste0("RMSE (", ts_unit[[sc]], ")"),
         PBIAS = "PBIAS (%)",
         KGE   = "KGE (-)")
}

source_cols <- c("CHIRPS v2"       = "darkred",
                 "CHIRPS v3 (RNL)" = "darkblue",
                 "CHIRPS v3 (SAT)" = "darkgreen")

n_col <- length(levels(metrics$Metric))
n_row <- length(timescales)

plots <- list()
k <- 0
for (sc in timescales) {            # rows
  for (m in levels(metrics$Metric)) {   # columns
    k <- k + 1
    # Repeated-measures design: keep only stations with a valid value for ALL
    # products, and sort by Source then Station so that rows are aligned
    # across groups (required by ggwithinstats).
    d <- metrics |>
      filter(Timescale == sc, Metric == m) |>
      group_by(Station) |>
      filter(sum(is.finite(value)) == length(products)) |>
      ungroup() |>
      arrange(Source, Station) |>
      droplevels()
    
    ttl <- paste0(letters[k], ") ", metric_names[[m]], " - ", sc)
    
    # x-axis text (product names) only on the bottom-most row
    hide_x <- if (k > (n_row - 1) * n_col) theme() else
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    plots[[k]] <- ggwithinstats(
      data = d, x = Source, y = value,
      type = "robust",                 # trimmed means (robust repeated measures)
      pairwise.display = "none",
      bf.message = FALSE,
      xlab = NULL,
      ylab = ylab_fun(m, sc),
      title = ttl
    ) +
      scale_color_manual(values = source_cols) +
      scale_fill_manual(values = source_cols) +
      theme(plot.title    = element_text(colour = "black", size = 16,
                                         face = "bold", hjust = 0),
            plot.subtitle = element_text(colour = "black", size = 12,
                                         face = "italic", hjust = 0),
            axis.text     = element_text(colour = "black", size = 12),
            axis.title    = element_text(colour = "black", size = 12),
            axis.title.x  = element_blank()) +   # remove the "Source" title
      hide_x
  }
}

fig <- ggarrange(plotlist = plots, ncol = n_col, nrow = n_row)

# The statistical subtitle (about 8 in long at 12 pt) must fit in one panel, so
# each panel needs to be about 9.5 in wide. The canvas is sized for that, with no
# `scale` argument (scale < 1 shrinks the canvas and crops the subtitles).
# If a subtitle is still cropped, increase fig_width.
fig_width  <- 38
fig_height <- 26
ggsave("figures/2_chirps_metrics_multiscale.png", fig,
       width = fig_width, height = fig_height, units = "in", dpi = 500, bg = "white")