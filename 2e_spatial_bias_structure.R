root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# Spatial structure of the mean annual rainfall bias (1998-2020) of CHIRPS v2,
# CHIRPS v3 (RNL) and CHIRPS v3 (SAT).
#   1. Bias rasters (product - reference) produced by 2d_annual_bias_maps.R
#   2. Local Moran's I (LISA) of the bias for each product
#      (distance-based weights, permutation test, FDR-adjusted p-values)
#   3. Empirical semivariogram of the bias + fitted model for each product
#   4. Figure: 2 rows x 2 columns
#        a) CHIRPS v2 LISA map   b) CHIRPS v3 (RNL) LISA map
#        c) CHIRPS v3 (SAT) LISA map   d) semivariograms (3 products overlaid)
#      -> "figures/2_chirps_bias_spatial_structure.png"
#   5. Table for the paper -> "tables/2_chirps_bias_spatial_structure.csv"
#
# NB: "High"/"Low" in LISA refer to the bias relative to the mean bias of the
#     product (not relative to zero).
# =============================================================================

# ---- 0. Packages -------------------------------------------------------------
pkgs <- c("terra", "sf", "spdep", "gstat", "dplyr", "ggplot2", "ggspatial",
          "ggpubr", "scales")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 1. Settings -------------------------------------------------------------
outline_file <- "gis/shp/bf_deg.shp"
derived_dir  <- "data/derived"

products <- c("CHIRPS v2"       = "chirps_v2",
              "CHIRPS v3 (RNL)" = "chirps_v3_rnl",
              "CHIRPS v3 (SAT)" = "chirps_v3_sat")

# Colours of the products (figure conventions)
prod_cols <- c("CHIRPS v2"       = "darkred",
               "CHIRPS v3 (RNL)" = "darkblue",
               "CHIRPS v3 (SAT)" = "darkgreen")
prod_shapes <- c("CHIRPS v2" = 16, "CHIRPS v3 (RNL)" = 17, "CHIRPS v3 (SAT)" = 15)
prod_ltys   <- c("CHIRPS v2" = "solid", "CHIRPS v3 (RNL)" = "solid",
                 "CHIRPS v3 (SAT)" = "dashed")

# Projected CRS used to compute distances in km (UTM 30N)
crs_proj <- 32630

# LISA settings
lisa_dist_km  <- 30        # neighbourhood radius (km); cells are ~5 km wide
lisa_nsim     <- 499       # number of permutations
lisa_alpha    <- 0.05      # significance level (after p-value adjustment)
lisa_p_adjust <- "fdr"     # multiple-testing correction ("none" to disable)

# Semivariogram settings
vgm_cutoff_km <- 300                     # maximum distance (km)
vgm_width_km  <- 10                      # lag width (km)
vgm_trend     <- FALSE                   # TRUE = remove a latitudinal trend first
vgm_models    <- c("Sph", "Exp", "Gau")  # candidate models (best SSErr kept)

lisa_levels <- c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant")
lisa_cols   <- c("High-High"       = "#2166ac",   # same palette as the bias maps
                 "Low-Low"         = "#b35806",
                 "High-Low"        = "#92c5de",
                 "Low-High"        = "#fdb863",
                 "Not significant" = "grey90")

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Study area and bias rasters -----------------------------------------
outline <- st_read(outline_file, quiet = TRUE) |> st_transform(4326)

bb <- st_bbox(outline)
xlim_map <- c(bb[["xmin"]] - 0.05, bb[["xmax"]] + 0.05)
ylim_map <- c(bb[["ymin"]] - 0.05, bb[["ymax"]] + 0.05)

bias_files <- file.path(derived_dir, paste0("2_bias_annual_rainfall_", products, ".tif"))
if (!all(file.exists(bias_files))) {
  stop("Bias rasters not found in '", derived_dir, "'. Run 2d_annual_bias_maps.R first.")
}
bias <- setNames(lapply(bias_files, rast), names(products))

# Bias values (one row per cell) with projected coordinates in km
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

# ---- 3. Local Moran's I (LISA) ----------------------------------------------
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

# ---- 4. Semivariograms -------------------------------------------------------
fm_v <- if (vgm_trend) bias ~ Y else bias ~ 1

fit_vgm <- function(v_emp) {
  fit_one <- function(m) {
    init <- gstat::vgm(psill  = max(max(v_emp$gamma) - min(v_emp$gamma), 1e-6),
                       model  = m,
                       range  = max(v_emp$dist) / 2,
                       nugget = min(v_emp$gamma))
    tryCatch(suppressWarnings(gstat::fit.variogram(v_emp, init)),
             error = function(e) NULL)
  }
  fits <- Filter(Negate(is.null), setNames(lapply(vgm_models, fit_one), vgm_models))
  if (length(fits) == 0) stop("Variogram fitting failed for all candidate models.")
  sse <- sapply(fits, function(f) attr(f, "SSErr"))
  fits[[names(which.min(sse))]]
}

vg <- lapply(names(products), function(pr) {
  message("Semivariogram: ", pr, " ...")
  v_emp <- gstat::variogram(fm_v, locations = ~ X + Y, data = cells[[pr]],
                            cutoff = vgm_cutoff_km, width = vgm_width_km)
  f     <- fit_vgm(v_emp)
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

# ---- 5. Summary table --------------------------------------------------------
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

# ---- 5b. Table for the paper (tables/) ---------------------------------------
# One row per product: bias magnitude, global spatial autocorrelation (Moran's I),
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

# ---- 6. Figure ---------------------------------------------------------------
theme_panel <- theme_bw() +
  theme(plot.title   = element_text(colour = "black", size = 18,
                                    face = "bold", hjust = 0),
        axis.text    = element_text(colour = "black", size = 12),
        axis.title   = element_text(colour = "black", size = 12),
        legend.text  = element_text(colour = "black", size = 12),
        legend.title = element_text(colour = "black", size = 12),
        panel.grid.minor = element_blank())

# --- (a-c) LISA maps
map_plots <- list()
for (k in seq_along(products)) {
  pr <- names(products)[k]
  d  <- lisa[[pr]]$data
  map_plots[[k]] <- ggplot() +
    geom_raster(data = d, aes(x = x, y = y, fill = cluster)) +
    geom_sf(data = outline, fill = NA, colour = "black", linewidth = 0.4) +
    scale_fill_manual(values = lisa_cols, limits = lisa_levels, drop = FALSE,
                      name = "LISA cluster",
                      guide = guide_legend(nrow = 1, title.position = "left",
                                           title.vjust = 0.5)) +
    coord_sf(xlim = xlim_map, ylim = ylim_map, expand = FALSE) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           style = north_arrow_minimal(text_size = 10),
                           height = unit(1.2, "cm"), width = unit(1.2, "cm")) +
    annotation_scale(location = "br", unit_category = "metric",
                     width_hint = 0.3, text_cex = 0.9) +
    labs(title = paste0(letters[k], ") ", pr, " - LISA clusters"),
         x = "Longitude (°)", y = "Latitude (°)") +
    theme_panel +
    theme(legend.position = "bottom")
}

# Common legend taken from the first map
leg <- as_ggplot(get_legend(map_plots[[1]]))
map_plots <- lapply(map_plots, function(p) p + theme(legend.position = "none"))

# x-axis titles only on the bottom panels, y-axis title only in the left column
map_plots[[1]] <- map_plots[[1]] + theme(axis.title.x = element_blank())
map_plots[[2]] <- map_plots[[2]] + theme(axis.title.x = element_blank(),
                                         axis.title.y = element_blank())

# --- (d) semivariograms of the three products overlaid
lab_prod <- setNames(
  sprintf("%s (range = %.0f km)", names(products),
          sapply(vg, function(v) v$eff_range)),
  names(products))

p_d <- ggplot() +
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
  labs(title = "d) Bias semivariogram",
       x = "Distance (km)",
       y = expression(Semivariance~(mm^2~year^{-2}))) +
  theme_panel +
  theme(legend.position = "inside",
        legend.position.inside = c(0.98, 0.04),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black", fill = alpha("white", 0.85)),
        legend.key.width = unit(1.2, "cm"))

# --- Assembly: 2 rows x 2 columns + common LISA legend at the bottom
body <- ggarrange(map_plots[[1]], map_plots[[2]], map_plots[[3]], p_d,
                  ncol = 2, nrow = 2, widths = c(1, 1))
fig  <- ggarrange(body, leg, ncol = 1, nrow = 2, heights = c(1, 0.06))

ggsave("figures/2_chirps_bias_spatial_structure.png", fig,
       width = 16, height = 13, units = "in", dpi = 350, bg = "white",
       scale = 0.6)