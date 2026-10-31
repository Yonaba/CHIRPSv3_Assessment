root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# Bias of mean annual rainfall (1998-2020) of CHIRPS v2, CHIRPS v3 (RNL) and
# CHIRPS v3 (SAT) against a reference raster interpolated from stations.
#   1. Mean annual rainfall rasters of the three products (v3 aligned on v2)
#   2. Mean annual rainfall per station (observations) -> IDW interpolation on
#      the CHIRPS v2 grid (kriging, gstat) -> reference raster
#   3. Bias = product - reference (mm/year)
#   4. Figure: 3 rows (products) x 2 columns (bias map, bias density)
#      -> "figures/2_chirps_annual_bias_maps.png"
# =============================================================================

# ---- 0. Packages -------------------------------------------------------------
pkgs <- c("terra", "sf", "gstat", "dplyr", "tidyr", "lubridate", "ggplot2",
          "ggspatial", "ggpubr", "readr", "scales")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 1. Settings -------------------------------------------------------------
start_date <- as.Date("1998-01-01")
end_date   <- as.Date("2020-12-31")
years      <- year(start_date):year(end_date)

stations_file <- "data/csv/p_stations.csv"
obs_file      <- "data/csv/p_stations_daily_data.csv"
chirps_dir    <- "data/chirps"
outline_file  <- "gis/shp/bf_deg.shp"

products <- c("CHIRPS v2"       = "chirps_v2",
              "CHIRPS v3 (RNL)" = "chirps_v3_rnl",
              "CHIRPS v3 (SAT)" = "chirps_v3_sat")

# Station annual totals: a year is kept if at least this fraction of its days
# is observed (annual total is then scaled to the full year); a station is kept
# if it has at least `min_years` valid years.
min_valid_frac <- 1.0
min_years      <- 10

# Kriging settings
krig_trend <- FALSE                          # TRUE = universal kriging (latitude trend)
vgm_models <- c("Sph", "Exp", "Gau")         # candidate variogram models (best SSErr kept)

dir.create("figures", showWarnings = FALSE)
dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)

# ---- 2. Study area -----------------------------------------------------------
outline <- st_read(outline_file, quiet = TRUE) |> st_transform(4326)
outline_v <- vect(outline)

bb <- st_bbox(outline)
xlim_map <- c(bb[["xmin"]] - 0.05, bb[["xmax"]] + 0.05)
ylim_map <- c(bb[["ymin"]] - 0.05, bb[["ymax"]] + 0.05)

# Extent used to crop the rasters (outline extent + margin for resampling)
crop_ext <- ext(bb[["xmin"]] - 0.25, bb[["xmax"]] + 0.25,
                bb[["ymin"]] - 0.25, bb[["ymax"]] + 0.25)

# ---- 3. Mean annual rainfall of the CHIRPS products -------------------------
mean_annual_chirps <- function(folder) {
  files <- list.files(file.path(chirps_dir, folder), pattern = "[0-9]{4}\\.tif$",
                      full.names = TRUE)
  yrs <- as.integer(sub(".*([0-9]{4})\\.tif$", "\\1", files))
  files <- files[yrs %in% years]; yrs <- yrs[yrs %in% years]
  if (length(files) == 0) stop("No tif files found in ", folder)
  if (!all(years %in% yrs)) warning("Missing years in ", folder, ": ",
                                    paste(setdiff(years, yrs), collapse = ", "))
  
  annual <- lapply(files, function(f) {
    r <- crop(rast(f), crop_ext, snap = "out")
    r[r < 0] <- NA                              # NoData (-9999)
    app(r, "sum", na.rm = FALSE)                # annual total (mm/year)
  })
  app(rast(annual), "mean")                     # mean over years
}

mean_rain <- lapply(products, function(f) {
  message("Mean annual rainfall: ", f, " ...")
  mean_annual_chirps(f)
})

# Template = CHIRPS v2 grid cropped to the study area; v3 aligned on it
template <- crop(mean_rain[["CHIRPS v2"]], ext(outline_v), snap = "out")
mean_rain[["CHIRPS v2"]] <- template
for (pr in setdiff(names(products), "CHIRPS v2")) {
  mean_rain[[pr]] <- resample(mean_rain[[pr]], template, method = "bilinear")
}

# Keep only the study area
mean_rain <- lapply(mean_rain, function(r) mask(r, outline_v))

# ---- 4. Mean annual rainfall at stations and interpolated reference ---------
stations <- read_csv(stations_file, show_col_types = FALSE) |>
  mutate(Station = trimws(as.character(Station))) |>
  filter(!is.na(Latitude), !is.na(Longitude))

obs_wide <- read_csv(obs_file, show_col_types = FALSE, name_repair = "minimal")
names(obs_wide) <- trimws(names(obs_wide))
names(obs_wide)[1] <- "Date"

obs_long <- obs_wide |>
  mutate(Date = as.Date(as.character(Date), format = "%m/%d/%Y")) |>
  filter(Date >= start_date, Date <= end_date) |>
  mutate(across(-Date, as.numeric)) |>
  pivot_longer(-Date, names_to = "Station", values_to = "obs") |>
  filter(Station %in% stations$Station)

st_annual <- obs_long |>
  mutate(Year = year(Date),
         ndays = ifelse(leap_year(Year), 366, 365)) |>
  group_by(Station, Year, ndays) |>
  summarise(n_valid = sum(is.finite(obs)),
            total = sum(obs, na.rm = TRUE), .groups = "drop") |>
  filter(n_valid / ndays >= min_valid_frac) |>
  mutate(total = total * ndays / n_valid)       # scale to full year if gaps

st_mean <- st_annual |>
  group_by(Station) |>
  summarise(n_years = n(), P = mean(total), .groups = "drop") |>
  filter(n_years >= min_years) |>
  left_join(stations |> select(Station, Longitude, Latitude), by = "Station")

message(nrow(st_mean), " stations used for the reference raster.")

st_pts <- st_as_sf(st_mean, coords = c("Longitude", "Latitude"),
                   crs = 4326, remove = FALSE)

# Kriging model: ordinary (P ~ 1) or universal with a latitude trend
fm <- if (krig_trend) P ~ Latitude else P ~ 1

# Empirical variogram (distances in km for lon/lat data) and model selection
v_emp <- gstat::variogram(fm, st_pts)
fit_one <- function(m) {
  init <- gstat::vgm(psill   = max(max(v_emp$gamma) - min(v_emp$gamma), 1e-6),
                     model   = m,
                     range   = max(v_emp$dist) / 2,
                     nugget  = min(v_emp$gamma))
  tryCatch(suppressWarnings(gstat::fit.variogram(v_emp, init)),
           error = function(e) NULL)
}
fits <- Filter(Negate(is.null), setNames(lapply(vgm_models, fit_one), vgm_models))
if (length(fits) == 0) stop("Variogram fitting failed for all candidate models.")
sse <- sapply(fits, function(f) attr(f, "SSErr"))
best <- names(which.min(sse))
vgm_fit <- fits[[best]]
message("Selected variogram model: ", best)
print(vgm_fit)

# Prediction on the cells of the CHIRPS v2 grid located inside the study area
cell_ids <- which(!is.na(values(mask(template, outline_v), mat = FALSE)))
xy <- xyFromCell(template, cell_ids)
newdata <- st_as_sf(data.frame(Longitude = xy[, 1], Latitude = xy[, 2]),
                    coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

krig <- gstat::krige(fm, st_pts, newdata, model = vgm_fit)

reference <- template; values(reference) <- NA
reference[cell_ids] <- pmax(krig$var1.pred, 0)     # rainfall cannot be negative
reference <- mask(reference, outline_v)
names(reference) <- "reference"

krig_var <- template; values(krig_var) <- NA
krig_var[cell_ids] <- krig$var1.var
krig_var <- mask(krig_var, outline_v)
names(krig_var) <- "kriging_variance"
writeRaster(krig_var, "data/derived/2_reference_kriging_variance.tif", overwrite = TRUE)

# ---- 5. Bias rasters (product - reference, mm/year) -------------------------
bias <- lapply(mean_rain, function(r) {
  b <- r - reference
  names(b) <- "bias"
  b
})

# Save derived rasters
writeRaster(reference, "data/derived/2_reference_annual_rainfall.tif", overwrite = TRUE)
for (pr in names(products)) {
  writeRaster(mean_rain[[pr]],
              paste0("data/derived/2_mean_annual_rainfall_", products[[pr]], ".tif"),
              overwrite = TRUE)
  writeRaster(bias[[pr]],
              paste0("data/derived/2_bias_annual_rainfall_", products[[pr]], ".tif"),
              overwrite = TRUE)
}

# ---- 6. Figure ---------------------------------------------------------------
bias_df <- lapply(bias, function(b) as.data.frame(b, xy = TRUE, na.rm = TRUE))

all_bias <- unlist(lapply(bias_df, function(d) d$bias))
lim      <- as.numeric(quantile(abs(all_bias), 0.99))     # symmetric map range
x_dens   <- as.numeric(quantile(all_bias, c(0.005, 0.995))) # common density x-range

theme_panel <- theme_bw() +
  theme(plot.title   = element_text(colour = "black", size = 18,
                                    face = "bold", hjust = 0),
        axis.text    = element_text(colour = "black", size = 12),
        axis.title   = element_text(colour = "black", size = 12),
        legend.text  = element_text(colour = "black", size = 12),
        legend.title = element_text(colour = "black", size = 12),
        panel.grid.minor = element_blank())

plots <- list()
k <- 0
for (pr in names(products)) {
  d <- bias_df[[pr]]
  
  # --- (left) bias map
  k <- k + 1
  plots[[k]] <- ggplot() +
    geom_raster(data = d, aes(x = x, y = y, fill = bias)) +
    geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.4) +
    scale_fill_gradient2(
      low = "#b35806", mid = "white", high = "#2166ac", midpoint = 0,
      limits = c(-lim, lim), oob = scales::squish,
      name = "Bias (mm/year)",
      guide = guide_colourbar(direction = "horizontal",
                              title.position = "left", title.vjust = 0.8,
                              barwidth = unit(10, "cm"),
                              barheight = unit(0.6, "cm"))) +
    coord_sf(xlim = xlim_map, ylim = ylim_map, expand = FALSE) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           style = north_arrow_minimal(text_size = 10),
                           height = unit(1.2, "cm"), width = unit(1.2, "cm")) +
    annotation_scale(location = "br", unit_category = "metric",
                     width_hint = 0.3, text_cex = 0.9) +
    labs(title = paste0(letters[k], ") ", pr, " - Bias map"),
         x = "Longitude (°)", y = "Latitude (°)") +
    theme_panel +
    theme(legend.position = "bottom")
  
  # --- (right) bias density
  k <- k + 1
  med <- median(d$bias)
  plots[[k]] <- ggplot(d, aes(x = bias)) +
    geom_density(fill = "grey70", colour = "black", alpha = 0.7, linewidth = 0.5) +
    geom_vline(xintercept = med, colour = "darkred", linetype = "dashed",
               linewidth = 0.8) +
    annotate("text", x = med, y = Inf,
             label = sprintf("Median = %.1f", med),
             colour = "darkred", size = 5, hjust = -0.05, vjust = 1.5) +
    coord_cartesian(xlim = x_dens) +
    labs(title = paste0(letters[k], ") "," Bias distribution"),
         x = "Bias (mm/year)", y = "Density") +
    theme_panel +
    theme(legend.position = "none")
  
  # x-axis titles only on the bottom-most panels (both columns)
  if (pr != tail(names(products), 1)) {
    plots[[k - 1]] <- plots[[k - 1]] + theme(axis.title.x = element_blank())
    plots[[k]]     <- plots[[k]]     + theme(axis.title.x = element_blank())
  }
}

# First column (maps) 1.5 times wider than the second (density plots)
fig <- ggarrange(plotlist = plots, ncol = 2, nrow = 3, widths = c(1.5, 1.2),
                 common.legend = TRUE, legend = "bottom")

ggsave("figures/2_chirps_annual_bias_maps.png", fig,
       width = 16, height = 20, units = "in", dpi = 350, bg = "white",
       scale = 0.6)
