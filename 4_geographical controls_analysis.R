root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# Geographical controls analysis of CHIRPS v2, CHIRPS v3 (RNL) and CHIRPS v3 (SAT)
# -----------------------------------------------------------------------------
# Question: do altitude and/or latitude affect the accuracy of the products?
#
# Method
#   1. For each station, extract the daily CHIRPS series of the pixel containing
#      the station (3 products) and pair it with observed daily rainfall.
#   2. Per station and product, compute the relative bias (%) on paired valid
#      days of the 1998-2020 period:
#        100 * (sum(CHIRPS) - sum(obs)) / sum(obs)
#   3. Altitude = mean of the HydroSHEDS DEM over the CHIRPS pixel footprint
#      (consistent with the CHIRPS resolution); latitude from the station file.
#   4. Figure: 2 x 3 panels (rows: altitude, latitude; columns: products) with
#      linear fit + 95% CI, R2, Spearman rho, slope and p-value.
#   5. Table: simple and multiple regressions (bias ~ altitude + latitude) to
#      isolate confounding between altitude and latitude.
#
# Outputs
#   figures/4_orography_relative_bias.png
#   tables/3_multiple_regression_bias.csv
#   tables/3_station_bias_altitude_latitude.csv   (station-level inputs)
#
# The script sets the working directory to its own folder, which must be the
# project root (the folder containing data/, gis/, figures/, tables/).
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ---- 0. Settings -------------------------------------------------------------
f_stations <- "data/csv/p_stations.csv"
f_obs      <- "data/csv/p_stations_daily_data.csv"
dir_chirps <- "data/chirps"
f_dem      <- "gis/raster/hydrosheds_dem.tif"
dir_fig    <- "figures"
dir_tab    <- "tables"

year_start <- 1998
year_end   <- 2020
min_cov    <- 0.80  # min. fraction of valid paired days for a station-year to be kept
min_years  <- 10    # min. number of kept years for a station to be analysed

products <- data.frame(
  id     = c("chirps_v2", "chirps_v3_rnl", "chirps_v3_sat"),
  label  = c("CHIRPS v2", "CHIRPS v3 (RNL)", "CHIRPS v3 (SAT)"),
  colour = c("#8B0000", "#00008B", "#006400"),  # dark red, dark blue, dark green
  stringsAsFactors = FALSE
)

dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)

days_in_year <- function(y) ifelse((y %% 4 == 0 & y %% 100 != 0) | y %% 400 == 0, 366, 365)

# ---- 1. Observations (wide -> long) -----------------------------------------
stations <- read.csv(f_stations, check.names = FALSE, strip.white = TRUE,
                     stringsAsFactors = FALSE)
stations$Station <- trimws(stations$Station)

obs_wide <- read.csv(f_obs, check.names = FALSE, strip.white = TRUE,
                     stringsAsFactors = FALSE)
names(obs_wide) <- trimws(names(obs_wide))
names(obs_wide)[1] <- "Date"
obs_wide$Date <- as.Date(obs_wide$Date, format = "%m/%d/%Y")
yr <- as.integer(format(obs_wide$Date, "%Y"))
obs_wide <- obs_wide[!is.na(yr) & yr >= year_start & yr <= year_end, ]

obs <- obs_wide |>
  pivot_longer(-Date, names_to = "Station", values_to = "obs") |>
  mutate(obs = suppressWarnings(as.numeric(obs)),
         obs = ifelse(obs < 0, NA, obs))  # negative values = missing-data codes

obs_stations <- setdiff(names(obs_wide), "Date")
if (length(setdiff(stations$Station, obs_stations)))
  message("Stations without daily data (dropped): ",
          paste(setdiff(stations$Station, obs_stations), collapse = ", "))
if (length(setdiff(obs_stations, stations$Station)))
  message("Daily-data columns without coordinates (dropped): ",
          paste(setdiff(obs_stations, stations$Station), collapse = ", "))

st <- stations |>
  filter(Station %in% obs_stations, !is.na(Latitude), !is.na(Longitude))

# ---- 2. Index the CHIRPS files ----------------------------------------------
index_files <- function(p) {
  f <- list.files(file.path(dir_chirps, p), pattern = "[0-9]{4}\\.tif$",
                  full.names = TRUE, ignore.case = TRUE)
  d <- data.frame(file = f,
                  year = as.integer(sub(".*([0-9]{4})\\.tif$", "\\1", f, ignore.case = TRUE)))
  d <- d[d$year >= year_start & d$year <= year_end, ]
  d[order(d$year), ]
}
chirps_files <- setNames(lapply(products$id, index_files), products$id)
for (p in products$id) {
  if (nrow(chirps_files[[p]]) == 0) stop("No tif file found for ", p)
  miss <- setdiff(year_start:year_end, chirps_files[[p]]$year)
  if (length(miss)) warning(p, ": missing years: ", paste(miss, collapse = ", "))
}

# ---- 3. Altitude: DEM mean over the CHIRPS pixel footprint -------------------
ref  <- rast(chirps_files[[1]]$file[1])[[1]]
for (p in products$id[-1]) {
  r1 <- rast(chirps_files[[p]]$file[1])[[1]]
  if (!isTRUE(compareGeom(ref, r1, stopOnError = FALSE)))
    warning(p, " is not on the same grid as ", products$id[1],
            "; the altitude of the first grid's pixel footprint is used for all products.")
}

xy    <- as.matrix(st[, c("Longitude", "Latitude")])
cells <- cellFromXY(ref, xy)
if (anyNA(cells)) {
  message("Stations outside the CHIRPS extent (dropped): ",
          paste(st$Station[is.na(cells)], collapse = ", "))
  st <- st[!is.na(cells), ]; xy <- xy[!is.na(cells), , drop = FALSE]; cells <- cells[!is.na(cells)]
}

cen   <- xyFromCell(ref, cells)
half  <- res(ref) / 2
polys <- do.call(rbind, lapply(seq_len(nrow(cen)), function(i)
  as.polygons(ext(cen[i, 1] - half[1], cen[i, 1] + half[1],
                  cen[i, 2] - half[2], cen[i, 2] + half[2]), crs = "EPSG:4326")))

dem <- rast(f_dem)  # EPSG:32630 (UTM 30N, metres)

# Station points and pixel footprints are in EPSG:4326: reproject them to the DEM CRS
polys_dem <- project(polys, crs(dem))
pts_dem   <- project(vect(xy, crs = "EPSG:4326"), crs(dem))

dem <- crop(dem, ext(polys_dem) + 2000)  # 2 km margin; speeds up the extraction
st$elev_pixel <- terra::extract(dem, polys_dem, fun = mean, na.rm = TRUE)[, 2]
st$elev_point <- terra::extract(dem, pts_dem)[, 2]
st <- st[!is.na(st$elev_pixel), ]
message(nrow(st), " stations retained for the analysis.")

# ---- 4. Extract daily CHIRPS at the stations ---------------------------------
pts_chirps <- vect(xy, crs = "EPSG:4326")  # station points (same order as st)

extract_chirps <- function(fdf) {
  bind_rows(lapply(seq_len(nrow(fdf)), function(i) {
    r <- rast(fdf$file[i])
    y <- fdf$year[i]
    n <- nlyr(r)
    if (n != days_in_year(y))
      warning(basename(fdf$file[i]), ": ", n, " bands (expected ", days_in_year(y), ")")
    # Dates from the band names (e.g. "19980101"); fallback: consecutive days from 1 Jan
    dates <- as.Date(names(r), format = "%Y%m%d")
    if (anyNA(dates) || length(dates) != n)
      dates <- seq(as.Date(sprintf("%d-01-01", y)), by = "day", length.out = n)
    # ID = FALSE: no ID column to drop (dropping column 1 blindly removed band 1)
    m <- as.matrix(terra::extract(r, pts_chirps, ID = FALSE))  # stations x days
    if (ncol(m) != n || nrow(m) != nrow(st))
      stop(basename(fdf$file[i]), ": extracted ", nrow(m), " x ", ncol(m),
           " values, expected ", nrow(st), " x ", n)
    m[m < 0] <- NA                                              # no-data codes (e.g. -9999)
    data.frame(Station = rep(st$Station, times = n),
               Date    = rep(dates, each = nrow(st)),
               chirps  = as.vector(m))
  }))
}

chirps_long <- lapply(products$id, function(p) {
  message("Extracting ", p, " ...")
  extract_chirps(chirps_files[[p]])
})
names(chirps_long) <- products$id

# ---- 5. Relative bias per station and product ---------------------------------
station_bias <- function(cl) {
  paired <- obs |>
    inner_join(cl, by = c("Station", "Date")) |>
    filter(!is.na(obs), !is.na(chirps)) |>
    mutate(year = as.integer(format(Date, "%Y")))
  
  keep <- paired |>
    count(Station, year) |>
    filter(n / days_in_year(year) >= min_cov) |>
    select(Station, year)
  
  paired |>
    inner_join(keep, by = c("Station", "year")) |>
    group_by(Station) |>
    summarise(n_years = n_distinct(year), n_days = n(),
              obs_total = sum(obs), chirps_total = sum(chirps), .groups = "drop") |>
    filter(n_years >= min_years, obs_total > 0) |>
    mutate(rel_bias = 100 * (chirps_total - obs_total) / obs_total)
}

dat <- bind_rows(lapply(seq_len(nrow(products)), function(i) {
  station_bias(chirps_long[[products$id[i]]]) |> mutate(product = products$label[i])
})) |>
  left_join(st[, c("Station", "Longitude", "Latitude", "elev_pixel", "elev_point")],
            by = "Station") |>
  rename(lat = Latitude, elev = elev_pixel) |>
  mutate(product = factor(product, levels = products$label))

# write.csv(
#   dat |> select(product, Station, Longitude, lat, elev, elev_point,
#                 n_years, n_days, obs_total, chirps_total, rel_bias) |>
#     rename(Latitude = lat, Altitude_pixel_mean_m = elev, Altitude_point_m = elev_point) |>
#     mutate(across(where(is.numeric), \(x) round(x, 4))),
#   file.path(dir_tab, "3_station_bias_altitude_latitude.csv"), row.names = FALSE
# )

# ---- 6. Figure: 2 x 3 panels --------------------------------------------------
rng  <- range(c(dat$rel_bias, 0), na.rm = TRUE)
span <- diff(rng)
ylim <- c(rng[1] - 0.05 * span, rng[2] + 0.40 * span)  # head-room for the annotation

make_panel <- function(d, xvar, xlab, col, title = NULL, show_ylab = FALSE,
                       slope_scale = 1, slope_unit = "") {
  fit  <- lm(reformulate(xvar, "rel_bias"), data = d)
  sm   <- summary(fit)
  rho  <- suppressWarnings(cor(d[[xvar]], d$rel_bias, method = "spearman"))
  pval <- sm$coefficients[2, 4]
  ptxt <- if (pval < 0.001) "p < 0.001" else sprintf("p = %.3f", pval)
  lab  <- sprintf("R² = %.2f  |  Spearman ρ = %.2f\nslope = %.2f %s\n%s  (n = %d)",
                  sm$r.squared, rho, unname(coef(fit)[2]) * slope_scale, slope_unit,
                  ptxt, nrow(d))
  
  ggplot(d, aes(x = .data[[xvar]], y = rel_bias)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = col, fill = col, alpha = 0.15, linewidth = 0.9) +
    geom_point(colour = col, alpha = 0.85, size = 2.2) +
    annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.04, vjust = 1.2,
             size = 3.1, lineheight = 0.95, colour = "black") +
    coord_cartesian(ylim = ylim) +
    labs(title = title, x = xlab, y = if (show_ylab) "Relative bias (%)" else NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", colour = col, hjust = 0, size = 18),
          axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          axis.line = element_line(colour = "grey30"),
          axis.ticks = element_line(colour = "grey30"))
}

panels <- list()
for (i in seq_len(nrow(products))) {
  d   <- dat |> filter(product == products$label[i])
  col <- products$colour[i]
  panels[[i]]     <- make_panel(d, "elev", "Altitude (m a.s.l.)", col,
                                title = products$label[i], show_ylab = (i == 1),
                                slope_scale = 100, slope_unit = "% per 100 m")
  panels[[i + 3]] <- make_panel(d, "lat", "Latitude (°N)", col,
                                show_ylab = (i == 1),
                                slope_scale = 1, slope_unit = "% per °")
}

fig <- wrap_plots(panels, ncol = 3, nrow = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
ggsave(file.path(dir_fig, "4_geo_control_relative_bias.png"), fig,
       width = 11, height = 7, dpi = 500, device = png_dev, bg = "white")

# ---- 7. Regression table (simple + multiple, to isolate confounding) ---------
term_lab <- c("(Intercept)" = "Intercept",
              "elev100"     = "Altitude (per 100 m)",
              "lat"         = "Latitude (per degree)")

reg_one <- function(d, label) {
  if (nrow(d) < 5) return(NULL)
  d <- d |> mutate(elev100 = elev / 100)
  specs  <- list("Altitude only"       = "elev100",
                 "Latitude only"       = "lat",
                 "Altitude + Latitude" = c("elev100", "lat"))
  cor_el <- cor(d$elev100, d$lat)
  
  bind_rows(lapply(names(specs), function(m) {
    vars <- specs[[m]]
    f    <- lm(reformulate(vars, "rel_bias"), data = d)
    s    <- summary(f)
    cf   <- s$coefficients
    ci   <- confint(f)
    z    <- d
    z[c("rel_bias", vars)] <- lapply(z[c("rel_bias", vars)], \(x) as.numeric(scale(x)))
    fz   <- lm(reformulate(vars, "rel_bias"), data = z)
    
    data.frame(
      product   = label,
      model     = m,
      term      = unname(term_lab[rownames(cf)]),
      estimate  = cf[, 1],
      std_error = cf[, 2],
      ci_low    = ci[, 1],
      ci_high   = ci[, 2],
      t_value   = cf[, 3],
      p_value   = cf[, 4],
      std_beta  = c(NA, unname(coef(fz)[vars])),
      r2        = s$r.squared,
      adj_r2    = s$adj.r.squared,
      n_stations = nobs(f),
      vif       = if (length(vars) == 2) 1 / (1 - cor_el^2) else NA_real_,
      cor_altitude_latitude = cor_el,
      row.names = NULL
    )
  }))
}

reg_tab <- bind_rows(lapply(products$label, function(l)
  reg_one(dat |> filter(product == l), l))) |>
  mutate(across(where(is.numeric), \(x) signif(x, 4)))

write.csv(reg_tab, file.path(dir_tab, "3_multiple_regression_bias.csv"), row.names = FALSE)

print(as.data.frame(reg_tab), row.names = FALSE)
message("Done: figures/4_geo_control_relative_bias.png, tables/3_*.csv")
