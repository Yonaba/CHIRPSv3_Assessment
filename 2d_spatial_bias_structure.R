root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# Bias of mean annual rainfall (1998-2020) of CHIRPS v2, CHIRPS v3 (RNL) and
# CHIRPS v3 (SAT) against a station-based reference, and spatial structure of
# that bias.
#   1. Mean annual rainfall rasters of the three products (v3 aligned on v2)
#   2. Mean annual rainfall per station (observations) -> ordinary kriging
#      (gstat) on the CHIRPS v2 grid -> reference raster
#   3. Bias = product - reference (mm/year)
#   4. Local Moran's I (LISA) of the bias for each product
#      (distance-based weights, permutation test, FDR-adjusted p-values)
#   5. Empirical semivariogram of the bias + fitted model for each product
#   6. Tables -> "data/derived/2_bias_spatial_structure_summary.csv"
#               "tables/2_chirps_bias_spatial_structure.csv"
#   7. Figure: 2 columns x 3 rows + 1 full-width row
#        column 1: bias maps        a) CHIRPS v2   b) v3 (RNL)   c) v3 (SAT)
#        column 2: LISA clusters    d) CHIRPS v2   e) v3 (RNL)   f) v3 (SAT)
#        last row: g) bias semivariograms (3 products overlaid)
#      -> "figures/2_chirps_annual_bias_spatial_structure.png"
#
# NB: "High"/"Low" in LISA refer to the bias relative to the mean bias of the
#     product (not relative to zero).
# =============================================================================

# ---- 0. Packages -------------------------------------------------------------
pkgs <- c("terra", "sf", "gstat", "spdep", "dplyr", "tidyr", "lubridate",
          "ggplot2", "ggspatial", "ggpubr", "readr", "scales")
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
derived_dir   <- "data/derived"
fig_file      <- "figures/2_chirps_annual_bias_spatial_structure.png"

products <- c("CHIRPS v2"       = "chirps_v2",
              "CHIRPS v3 (RNL)" = "chirps_v3_rnl",
              "CHIRPS v3 (SAT)" = "chirps_v3_sat")

# Colours of the products (figure conventions)
prod_cols   <- c("CHIRPS v2"       = "darkred",
                 "CHIRPS v3 (RNL)" = "darkblue",
                 "CHIRPS v3 (SAT)" = "darkgreen")
prod_shapes <- c("CHIRPS v2" = 16, "CHIRPS v3 (RNL)" = 17, "CHIRPS v3 (SAT)" = 15)
prod_ltys   <- c("CHIRPS v2" = "solid", "CHIRPS v3 (RNL)" = "solid",
                 "CHIRPS v3 (SAT)" = "dashed")

# Station annual totals: a year is kept if at least this fraction of its days
# is observed (annual total is then scaled to the full year); a station is kept
# if it has at least `min_years` valid years.
min_valid_frac <- 1.0
min_years      <- 10

# Kriging of the reference raster
krig_trend <- FALSE                          # TRUE = universal kriging (latitude trend)

# Candidate variogram models (best SSErr kept) - kriging and bias semivariograms
vgm_models <- c("Sph", "Exp", "Gau")

# Projected CRS used to compute distances in km (UTM 30N)
crs_proj <- 32630

# LISA settings
lisa_dist_km  <- 30        # neighbourhood radius (km); cells are ~5 km wide
lisa_nsim     <- 499       # number of permutations
lisa_alpha    <- 0.05      # significance level (after p-value adjustment)
lisa_p_adjust <- "fdr"     # multiple-testing correction ("none" to disable)

# Bias semivariogram settings
vgm_cutoff_km <- 300       # maximum distance (km)
vgm_width_km  <- 10        # lag width (km)
vgm_trend     <- FALSE     # TRUE = remove a latitudinal trend first

lisa_levels <- c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant")
lisa_cols   <- c("High-High"       = "#2166ac",   # same palette as the bias maps
                 "Low-Low"         = "#b35806",
                 "High-Low"        = "#92c5de",
                 "Low-High"        = "#fdb863",
                 "Not significant" = "grey90")

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Study area -----------------------------------------------------------
outline   <- st_read(outline_file, quiet = TRUE) |> st_transform(4326)
outline_v <- vect(outline)

bb <- st_bbox(outline)
xlim_map <- c(bb[["xmin"]] - 0.05, bb[["xmax"]] + 0.05)
ylim_map <- c(bb[["ymin"]] - 0.05, bb[["ymax"]] + 0.05)

# Extent used to crop the rasters (outline extent + margin for resampling)
crop_ext <- ext(bb[["xmin"]] - 0.25, bb[["xmax"]] + 0.25,
                bb[["ymin"]] - 0.25, bb[["ymax"]] + 0.25)

# ---- 3. Helper: best-fitting variogram model ---------------------------------
fit_best_vgm <- function(v_emp, models = vgm_models) {
  fit_one <- function(m) {
    init <- gstat::vgm(psill  = max(max(v_emp$gamma) - min(v_emp$gamma), 1e-6),
                       model  = m,
                       range  = max(v_emp$dist) / 2,
                       nugget = min(v_emp$gamma))
    tryCatch(suppressWarnings(gstat::fit.variogram(v_emp, init)),
             error = function(e) NULL)
  }
  fits <- Filter(Negate(is.null), setNames(lapply(models, fit_one), models))
  if (length(fits) == 0) stop("Variogram fitting failed for all candidate models.")
  sse <- sapply(fits, function(f) attr(f, "SSErr"))
  fits[[names(which.min(sse))]]
}

# ---- 4. Mean annual rainfall of the CHIRPS products -------------------------
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

# ---- 5. Mean annual rainfall at stations and kriged reference ---------------
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
  left_join(stations |> dplyr::select(Station, Longitude, Latitude), by = "Station")

message(nrow(st_mean), " stations used for the reference raster.")

st_pts <- st_as_sf(st_mean, coords = c("Longitude", "Latitude"),
                   crs = 4326, remove = FALSE)

# Kriging model: ordinary (P ~ 1) or universal with a latitude trend
fm <- if (krig_trend) P ~ Latitude else P ~ 1

# Empirical variogram (distances in km for lon/lat data) and model selection
v_emp_ref <- gstat::variogram(fm, st_pts)
vgm_fit   <- fit_best_vgm(v_emp_ref)
message("Selected variogram model (reference): ", as.character(vgm_fit$model[2]))
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

# ---- 6. Bias rasters (product - reference, mm/year) -------------------------
bias <- lapply(mean_rain, function(r) {
  b <- r - reference
  names(b) <- "bias"
  b
})

# Save derived rasters
writeRaster(reference, file.path(derived_dir, "2_reference_annual_rainfall.tif"),
            overwrite = TRUE)
writeRaster(krig_var, file.path(derived_dir, "2_reference_kriging_variance.tif"),
            overwrite = TRUE)
for (pr in names(products)) {
  writeRaster(mean_rain[[pr]],
              file.path(derived_dir, paste0("2_mean_annual_rainfall_", products[[pr]], ".tif")),
              overwrite = TRUE)
  writeRaster(bias[[pr]],
              file.path(derived_dir, paste0("2_bias_annual_rainfall_", products[[pr]], ".tif")),
              overwrite = TRUE)
}

# ---- 7. Bias values per cell (lon/lat + projected coordinates in km) --------
get_cells <- function(b) {
  d <- as.data.frame(b, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "bias"
  xy_km <- st_as_sf(d, coords = c("x", "y"), crs = 4326, remove = FALSE) |>
    st_transform(crs_proj) |>
    st_coordinates()
  d$X <- xy_km[, 1] / 1000
  d$Y <- xy_km[, 2] / 1000
  d
}
cells <- lapply(bias, get_cells)

# ---- 8. Local Moran's I (LISA) ----------------------------------------------
run_lisa <- function(d) {
  nb <- dnearneigh(as.matrix(d[, c("X", "Y")]), d1 = 0, d2 = lisa_dist_km)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # Global Moran's I (same weights)
  mt <- moran.test(d$bias, lw, zero.policy = TRUE)
  
  # Local Moran's I with conditional permutations
  res <- as.data.frame(localmoran_perm(d$bias, lw, nsim = lisa_nsim,
                                       zero.policy = TRUE, iseed = 42))
  pcol <- grep("Sim$", colnames(res), value = TRUE)
  if (length(pcol) == 0) pcol <- grep("^Pr", colnames(res), value = TRUE)
  p     <- res[[pcol[1]]]
  p_adj <- p.adjust(p, method = lisa_p_adjust)
  
  # Moran scatterplot quadrants (standardised bias vs its spatial lag)
  z     <- as.numeric(scale(d$bias))
  lag_z <- lag.listw(lw, z, zero.policy = TRUE)
  
  cluster <- ifelse(p_adj >= lisa_alpha, "Not significant",
                    ifelse(z > 0 & lag_z > 0, "High-High",
                           ifelse(z < 0 & lag_z < 0, "Low-Low",
                                  ifelse(z > 0 & lag_z < 0, "High-Low", "Low-High"))))
  d$cluster <- factor(cluster, levels = lisa_levels)
  
  list(data = d,
       moran_I = unname(mt$estimate["Moran I statistic"]),
       moran_p = mt$p.value,
       mean_nb = mean(card(nb)))
}

lisa <- lapply(names(products), function(pr) {
  message("LISA: ", pr, " ...")
  run_lisa(cells[[pr]])
})
names(lisa) <- names(products)

# ---- 9. Bias semivariograms --------------------------------------------------
fm_v <- if (vgm_trend) bias ~ Y else bias ~ 1

vg <- lapply(names(products), function(pr) {
  message("Semivariogram: ", pr, " ...")
  v_emp <- gstat::variogram(fm_v, locations = ~ X + Y, data = cells[[pr]],
                            cutoff = vgm_cutoff_km, width = vgm_width_km)
  f     <- fit_best_vgm(v_emp)
  mod   <- as.character(f$model[2])
  rg    <- f$range[2]
  eff   <- switch(mod, Sph = rg, Exp = 3 * rg, Gau = sqrt(3) * rg, rg)
  vl    <- gstat::variogramLine(f, maxdist = vgm_cutoff_km, n = 200)
  message("  selected model: ", mod, " | effective range = ", round(eff), " km")
  list(emp  = data.frame(Product = pr, dist = v_emp$dist, gamma = v_emp$gamma),
       line = data.frame(Product = pr, dist = vl$dist, gamma = vl$gamma),
       model = mod, nugget = f$psill[1], psill = f$psill[2],
       range = rg, eff_range = eff)
})
names(vg) <- names(products)

vg_emp  <- bind_rows(lapply(vg, `[[`, "emp"))
vg_line <- bind_rows(lapply(vg, `[[`, "line"))
vg_emp$Product  <- factor(vg_emp$Product,  levels = names(products))
vg_line$Product <- factor(vg_line$Product, levels = names(products))

# ---- 10. Tables --------------------------------------------------------------
summary_tab <- bind_rows(lapply(names(products), function(pr) {
  d  <- lisa[[pr]]$data
  pc <- prop.table(table(d$cluster)) * 100
  v  <- vg[[pr]]
  data.frame(Product = pr,
             n_cells = nrow(d),
             mean_neighbours = round(lisa[[pr]]$mean_nb, 1),
             moran_I = round(lisa[[pr]]$moran_I, 3),
             pct_HH = round(pc[["High-High"]], 1),
             pct_LL = round(pc[["Low-Low"]], 1),
             pct_HL = round(pc[["High-Low"]], 1),
             pct_LH = round(pc[["Low-High"]], 1),
             pct_NS = round(pc[["Not significant"]], 1),
             vgm_model = v$model,
             nugget = round(v$nugget, 1),
             partial_sill = round(v$psill, 1),
             nugget_to_sill = round(v$nugget / (v$nugget + v$psill), 3),
             range_km = round(v$range, 1),
             effective_range_km = round(v$eff_range, 1))
}))
print(summary_tab)
write.csv(summary_tab, file.path(derived_dir, "2_bias_spatial_structure_summary.csv"),
          row.names = FALSE)

# Table for the paper: bias magnitude, global spatial autocorrelation (Moran's I),
# share of the territory in each LISA class, and semivariogram parameters.
# Spatial dependence class from the nugget-to-sill ratio (Cambardella et al., 1994):
#   < 0.25 strong, 0.25-0.75 moderate, > 0.75 weak.
dep_class <- function(r) ifelse(r < 0.25, "Strong", ifelse(r <= 0.75, "Moderate", "Weak"))
fmt_p     <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

paper_tab <- bind_rows(lapply(names(products), function(pr) {
  d  <- lisa[[pr]]$data
  pc <- prop.table(table(d$cluster)) * 100
  v  <- vg[[pr]]
  ns <- v$nugget / (v$nugget + v$psill)
  data.frame(
    "Product"                       = pr,
    "Cells (n)"                     = nrow(d),
    "Mean bias (mm/year)"           = round(mean(d$bias), 1),
    "Median bias (mm/year)"         = round(median(d$bias), 1),
    "SD of bias (mm/year)"          = round(sd(d$bias), 1),
    "Moran's I"                     = round(lisa[[pr]]$moran_I, 3),
    "Moran's I p-value"             = fmt_p(lisa[[pr]]$moran_p),
    "High-High (%)"                 = round(pc[["High-High"]], 1),
    "Low-Low (%)"                   = round(pc[["Low-Low"]], 1),
    "High-Low (%)"                  = round(pc[["High-Low"]], 1),
    "Low-High (%)"                  = round(pc[["Low-High"]], 1),
    "Not significant (%)"           = round(pc[["Not significant"]], 1),
    "Variogram model"               = v$model,
    "Nugget (mm2/year2)"            = round(v$nugget, 1),
    "Partial sill (mm2/year2)"      = round(v$psill, 1),
    "Nugget-to-sill ratio"          = round(ns, 2),
    "Spatial dependence"            = dep_class(ns),
    "Effective range (km)"          = round(v$eff_range, 0),
    check.names = FALSE, stringsAsFactors = FALSE)
}))
print(paper_tab)
write.csv(paper_tab, "tables/2_chirps_bias_spatial_structure.csv", row.names = FALSE)

# ---- 11. Figure --------------------------------------------------------------
n_prod <- length(products)

theme_panel <- theme_bw() +
  theme(plot.title   = element_text(colour = "black", size = 18,
                                    face = "bold", hjust = 0),
        axis.text    = element_text(colour = "black", size = 12),
        axis.title   = element_text(colour = "black", size = 12),
        legend.text  = element_text(colour = "black", size = 12),
        legend.title = element_text(colour = "black", size = 12),
        panel.grid.minor = element_blank())

# Common map decorations (outline, extent, north arrow, scale bar, labels)
map_decor <- function(title) {
  list(geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.4),
       coord_sf(xlim = xlim_map, ylim = ylim_map, expand = FALSE),
       annotation_north_arrow(location = "tr", which_north = "true",
                              style = north_arrow_minimal(text_size = 10),
                              height = unit(1.2, "cm"), width = unit(1.2, "cm")),
       annotation_scale(location = "br", unit_category = "metric",
                        width_hint = 0.3, text_cex = 0.9),
       labs(title = title, x = "Longitude (°)", y = "Latitude (°)"),
       theme_panel,
       theme(legend.position = "bottom"))
}

# --- Column 1 (a-c): bias maps
bias_df  <- lapply(bias, function(b) as.data.frame(b, xy = TRUE, na.rm = TRUE))
all_bias <- unlist(lapply(bias_df, function(d) d$bias))
lim      <- as.numeric(quantile(abs(all_bias), 0.99))     # symmetric map range

bias_maps <- lapply(seq_len(n_prod), function(k) {
  pr <- names(products)[k]
  ggplot() +
    geom_raster(data = bias_df[[pr]], aes(x = x, y = y, fill = bias)) +
    scale_fill_gradient2(
      low = "#b35806", mid = "white", high = "#2166ac", midpoint = 0,
      limits = c(-lim, lim), oob = scales::squish,
      name = "Bias (mm/year)",
      guide = guide_colourbar(direction = "horizontal",
                              title.position = "left", title.vjust = 0.8,
                              barwidth = unit(7, "cm"),
                              barheight = unit(0.6, "cm"))) +
    map_decor(paste0(letters[k], ") ", pr, " - Bias map"))
})

# --- Column 2 (d-f): LISA cluster maps
lisa_maps <- lapply(seq_len(n_prod), function(k) {
  pr <- names(products)[k]
  ggplot() +
    geom_raster(data = lisa[[pr]]$data, aes(x = x, y = y, fill = cluster)) +
    scale_fill_manual(values = lisa_cols, limits = lisa_levels, drop = FALSE,
                      name = "LISA cluster",
                      guide = guide_legend(nrow = 2, title.position = "left",
                                           title.vjust = 0.5)) +
    map_decor(paste0(letters[n_prod + k], ") ", pr, " - LISA clusters"))
})

# Column legends (placed in a dedicated row below the maps)
leg_bias <- as_ggplot(get_legend(bias_maps[[1]]))
leg_lisa <- as_ggplot(get_legend(lisa_maps[[1]]))

# Remove legends; x-axis titles only on the bottom row, y-axis title only in
# the left column
for (k in seq_len(n_prod)) {
  bias_maps[[k]] <- bias_maps[[k]] + theme(legend.position = "none")
  lisa_maps[[k]] <- lisa_maps[[k]] + theme(legend.position = "none",
                                           axis.title.y = element_blank())
  if (k < n_prod) {
    bias_maps[[k]] <- bias_maps[[k]] + theme(axis.title.x = element_blank())
    lisa_maps[[k]] <- lisa_maps[[k]] + theme(axis.title.x = element_blank())
  }
}

# --- Last row (g): semivariograms of the three products overlaid
lab_prod <- setNames(
  sprintf("%s (range = %.0f km)", names(products),
          sapply(vg, function(v) v$eff_range)),
  names(products))

p_vgm <- ggplot() +
  geom_line(data = vg_line,
            aes(x = dist, y = gamma, colour = Product, linetype = Product),
            linewidth = 1) +
  geom_point(data = vg_emp,
             aes(x = dist, y = gamma, colour = Product, shape = Product),
             size = 2.5) +
  scale_colour_manual(name = NULL, values = prod_cols,   labels = lab_prod) +
  scale_shape_manual( name = NULL, values = prod_shapes, labels = lab_prod) +
  scale_linetype_manual(name = NULL, values = prod_ltys, labels = lab_prod) +
  scale_x_continuous(limits = c(0, vgm_cutoff_km),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(title = paste0(letters[2 * n_prod + 1], ") Bias semivariogram"),
       x = "Distance (km)",
       y = expression(Semivariance~(mm^2~year^{-2}))) +
  theme_panel +
  theme(legend.position = "inside",
        legend.position.inside = c(0.99, 0.04),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black", fill = alpha("white", 0.85)),
        legend.key.width = unit(1.2, "cm"))

# --- Assembly
# Maps filled row by row: (bias v2, LISA v2), (bias RNL, LISA RNL), (bias SAT, LISA SAT)
map_grid   <- ggarrange(plotlist = unlist(Map(list, bias_maps, lisa_maps),
                                          recursive = FALSE),
                        ncol = 2, nrow = n_prod)
legend_row <- ggarrange(leg_bias, leg_lisa, ncol = 2)
fig        <- ggarrange(map_grid, legend_row, p_vgm,
                        ncol = 1, nrow = 3, heights = c(3, 0.18, 1.1))

ggsave(fig_file, fig,
       width = 18, height = 27, units = "in", dpi = 500, bg = "white",
       scale = 0.6)
message("Figure saved: ", fig_file)