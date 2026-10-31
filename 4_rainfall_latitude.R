root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# Latitude-time (Hovmoller-type) diagrams of the mean daily rainfall over
# Burkina Faso, 1998-2020: observations vs CHIRPS v2, CHIRPS v3 (RNL), CHIRPS v3 (SAT).
#   1. CHIRPS: daily rasters -> longitudinal mean for each latitude row (inside
#      the country outline) -> multi-year mean for each day of the year
#   2. Observations: multi-year mean daily rainfall at each station -> IDW
#      interpolation on the CHIRPS v2 grid -> same longitudinal mean per row
#   3. Light circular moving-average smoothing along the time axis
#   4. Figure: 2 rows x 2 columns (a Observations, b v2, c v3 RNL, d v3 SAT),
#      common colour scale (mm/day)
#      -> "figures/2_rainfall_latitude_time.png"
#
# Assumptions: each yearly CHIRPS tif holds 365/366 DAILY layers in calendar
# order; 29 February is dropped so that every year has 365 days.
# =============================================================================

# ---- 0. Packages -------------------------------------------------------------
pkgs <- c("terra", "sf", "dplyr", "tidyr", "lubridate", "ggplot2", "ggpubr",
          "readr", "scales")
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

# Same station selection as in 2d_annual_bias_maps.R: a year is used if at least
# `min_valid_frac` of its days are observed; a station is kept if it has at least
# `min_years` valid years.
min_valid_frac <- 1.0
min_years      <- 10

# Observation interpolation (IDW) and temporal smoothing
crs_proj    <- 32630      # UTM 30N, distances in km
idw_power   <- 2
smooth_days <- 15         # centred moving average (odd number of days; 1 = none)

# Colour classes (mm/day). Adjust to the range printed in the console.
brks <- c(0, 0.5, 1:10)
# Or let the classes follow the data (0 to the rounded-up maximum of all panels)
brks_auto <- FALSE
brks_step <- 0.5

dir.create("figures", showWarnings = FALSE)
dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)

# ---- 2. Study area and grid --------------------------------------------------
outline   <- st_read(outline_file, quiet = TRUE) |> st_transform(4326)
outline_v <- vect(outline)

bb <- st_bbox(outline)
crop_ext <- ext(bb[["xmin"]] - 0.25, bb[["xmax"]] + 0.25,
                bb[["ymin"]] - 0.25, bb[["ymax"]] + 0.25)

list_tifs <- function(folder) {
  f <- list.files(file.path(chirps_dir, folder), pattern = "[0-9]{4}\\.tif$",
                  full.names = TRUE)
  y <- as.integer(sub(".*([0-9]{4})\\.tif$", "\\1", f))
  f[y %in% years]
}

# Template = CHIRPS v2 grid cropped to the study area (as in 2d_annual_bias_maps.R)
f1 <- list_tifs(products[["CHIRPS v2"]])[1]
template <- crop(crop(rast(f1)[[1]], crop_ext, snap = "out"),
                 ext(outline_v), snap = "out")
nlat <- nrow(template)
lat  <- yFromRow(template, seq_len(nlat))

zone     <- mask(init(template, 1), outline_v)     # cells inside the country
cells_in <- which(!is.na(values(zone, mat = FALSE)))

# Calendar of a non-leap year (day-of-year index 1..365)
md_ref <- format(seq(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day"), "%m-%d")
n_doy  <- length(md_ref)

# ---- 3. CHIRPS: longitudinal mean per latitude row, mean over years ----------
chirps_lat_doy <- function(folder, is_ref) {
  files <- list_tifs(folder)
  if (length(files) == 0) stop("No tif files found in ", folder)
  acc <- NULL
  for (f in files) {
    r <- crop(rast(f), crop_ext, snap = "out")
    r[r < 0] <- NA                                        # NoData (-9999)
    if (!nlyr(r) %in% c(365, 366))
      stop(basename(f), ": ", nlyr(r), " layers (daily 365/366 layers expected).")
    if (nlyr(r) == 366) r <- r[[-60]]                     # drop 29 February
    r <- if (is_ref) crop(r, ext(template)) else resample(r, template, method = "bilinear")
    r <- mask(r, outline_v)
    arr <- as.array(r)                                    # rows x cols x days
    m <- rowMeans(aperm(arr, c(1, 3, 2)), dims = 2, na.rm = TRUE)  # rows x days
    m[is.nan(m)] <- NA
    acc <- if (is.null(acc)) m else acc + m
  }
  acc / length(files)
}

mats <- list()
for (pr in names(products)) {
  message("CHIRPS latitude-time matrix: ", pr, " ...")
  mats[[pr]] <- chirps_lat_doy(products[[pr]], is_ref = (pr == "CHIRPS v2"))
}

# ---- 4. Observations: station climatology -> IDW -> longitudinal mean --------
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
  filter(Station %in% stations$Station) |>
  mutate(Year = year(Date))

valid_sy <- obs_long |>
  mutate(ndays = ifelse(leap_year(Year), 366, 365)) |>
  group_by(Station, Year, ndays) |>
  summarise(n_valid = sum(is.finite(obs)), .groups = "drop") |>
  filter(n_valid / ndays >= min_valid_frac)

keep_st <- valid_sy |> count(Station) |> filter(n >= min_years) |> pull(Station)
message(length(keep_st), " stations used for the observed latitude-time diagram.")

clim_st <- obs_long |>
  semi_join(valid_sy |> filter(Station %in% keep_st), by = c("Station", "Year")) |>
  mutate(md = format(Date, "%m-%d")) |>
  filter(md != "02-29") |>
  mutate(doy = match(md, md_ref)) |>
  group_by(Station, doy) |>
  summarise(P = mean(obs, na.rm = TRUE), .groups = "drop")

# Matrix stations x days (mm/day)
V <- matrix(NA_real_, nrow = length(keep_st), ncol = n_doy)
V[cbind(match(clim_st$Station, keep_st), clim_st$doy)] <- clim_st$P
V[is.nan(V)] <- NA

st_info <- stations[match(keep_st, stations$Station), ]
to_km <- function(lon, lat) {
  st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326) |>
    st_transform(crs_proj) |> st_coordinates() / 1000
}
xy_cell <- xyFromCell(template, cells_in)
cell_km <- to_km(xy_cell[, 1], xy_cell[, 2])
st_km   <- to_km(st_info$Longitude, st_info$Latitude)

# IDW weights (cells x stations); missing station-days are ignored day by day
D <- sqrt(outer(cell_km[, 1], st_km[, 1], "-")^2 + outer(cell_km[, 2], st_km[, 2], "-")^2)
W <- 1 / pmax(D, 1)^idw_power
A  <- (!is.na(V)) * 1
V0 <- V; V0[is.na(V0)] <- 0
pred <- (W %*% V0) / (W %*% A)                            # cells x days

# Longitudinal mean for each latitude row
rows <- rowFromCell(template, cells_in)
zs   <- rowsum(pred, rows) / as.vector(table(rows))
M_obs <- matrix(NA_real_, nrow = nlat, ncol = n_doy)
M_obs[as.integer(rownames(zs)), ] <- zs
mats <- c(list("Observations" = M_obs), mats)

# ---- 5. Temporal smoothing and long format ----------------------------------
smooth_circ <- function(x, k) {
  if (k <= 1) return(x)
  h  <- (k - 1) %/% 2
  xx <- c(tail(x, h), x, head(x, h))
  as.numeric(stats::filter(xx, rep(1 / k, k), sides = 2))[(h + 1):(h + length(x))]
}
mats <- lapply(mats, function(M) t(apply(M, 1, smooth_circ, k = smooth_days)))

panel_names <- c("Observations", names(products))
clim <- bind_rows(lapply(panel_names, function(pr) {
  data.frame(Product = pr,
             lat = rep(lat, times = n_doy),
             doy = rep(seq_len(n_doy), each = nlat),
             P   = as.vector(mats[[pr]]))
})) |> filter(!is.na(P))
clim$Product <- factor(clim$Product, levels = panel_names)

#write.csv(clim, "data/derived/2_rainfall_latitude_time_climatology.csv", row.names = FALSE)

# Peak of each diagram (useful to adjust the colour classes)
peaks <- clim |>
  group_by(Product) |>
  slice_max(P, n = 1, with_ties = FALSE) |>
  transmute(Product, peak_mm_day = round(P, 1), latitude = round(lat, 2),
            date = format(as.Date(doy - 1, origin = "2001-01-01"), "%d %b"))
print(peaks)

# ---- 5b. Diagnostics of the rainfall range (mm/day) --------------------------
# Distribution of the plotted values, per panel
diag_tab <- clim |>
  group_by(Product) |>
  summarise(min = min(P), q25 = quantile(P, 0.25), median = median(P),
            q75 = quantile(P, 0.75), q95 = quantile(P, 0.95),
            q99 = quantile(P, 0.99), max = max(P), .groups = "drop") |>
  mutate(across(where(is.numeric), ~ round(.x, 2)))
print(as.data.frame(diag_tab))

# Share of plotted values falling in each colour class (are the top classes used?)
print(round(prop.table(table(clim$Product,
                             cut(clim$P, breaks = c(brks, Inf), right = FALSE,
                                 include.lowest = TRUE)), 1) * 100, 1))

# Quick look at the histograms
print(ggplot(clim, aes(x = P)) +
        geom_histogram(bins = 60, fill = "grey70", colour = "black") +
        facet_wrap(~ Product, ncol = 2) +
        labs(x = "Mean daily rainfall (mm/day)", y = "Count") +
        theme_bw())

# Colour classes following the data, if requested in the settings
if (brks_auto) brks <- seq(0, ceiling(max(clim$P) / brks_step) * brks_step, by = brks_step)

# ---- 6. Figure ---------------------------------------------------------------
n_b   <- length(brks)
lab_b <- c(paste0("<", brks[2]),
           paste(brks[2:(n_b - 1)], brks[3:n_b], sep = "-"),
           paste0(">", brks[n_b]))
pal   <- colorRampPalette(c("white", "#dfe7f5", "#a9b9e3", "#6a68b8", "#3f4fa3",
                            "#3c9a8f", "#4cb848", "#9bd23a", "#e6dc2e",
                            "#f3b02a", "#e8742a", "#d7301f"))(n_b)
clim$class <- cut(clim$P, breaks = c(brks, Inf), labels = lab_b,
                  include.lowest = TRUE, right = FALSE)

theme_panel <- theme_bw() +
  theme(plot.title   = element_text(colour = "black", size = 18,
                                    face = "bold", hjust = 0),
        axis.text    = element_text(colour = "black", size = 12),
        axis.title   = element_text(colour = "black", size = 12),
        legend.text  = element_text(colour = "black", size = 12),
        legend.title = element_text(colour = "black", size = 12),
        panel.grid.minor = element_blank())

month_start <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
lat_rng     <- range(lat[rowSums(!is.na(mats[[1]])) > 0])

plots <- list()
for (k in seq_along(panel_names)) {
  pr <- panel_names[k]
  plots[[k]] <- ggplot(filter(clim, Product == pr), aes(x = doy, y = lat, fill = class)) +
    geom_raster() +
    scale_fill_manual(values = setNames(pal, lab_b), drop = FALSE,
                      name = "Rainfall (mm/day)",
                      guide = guide_legend(nrow = 1, title.position = "left",
                                           title.vjust = 0.5, label.position = "bottom",
                                           keywidth = unit(1.3, "cm"),
                                           keyheight = unit(0.5, "cm"))) +
    scale_x_continuous(breaks = month_start, labels = paste("1", month.abb),
                       limits = c(1, n_doy), expand = c(0, 0)) +
    scale_y_continuous(labels = function(x) paste0(x, "°N"),
                       limits = lat_rng, expand = c(0, 0)) +
    labs(title = paste0(letters[k], ") ", pr), x = NULL, y = "Latitude") +
    theme_panel +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "bottom",
          legend.spacing.x = unit(0, "cm"))
  if (k %% 2 == 0) plots[[k]] <- plots[[k]] + theme(axis.title.y = element_blank())
  # x-axis text only on the bottom-most panels (here: panels c and d)
  if (k <= 2) plots[[k]] <- plots[[k]] + theme(axis.text.x = element_blank())
}

fig <- ggarrange(plotlist = plots, ncol = 2, nrow = 2,
                 common.legend = TRUE, legend = "bottom")

ggsave("figures/4_rainfall_latitude_time.png", fig,
       width = 24, height = 12, units = "in", dpi = 350, bg = "white",
       scale = 0.6)