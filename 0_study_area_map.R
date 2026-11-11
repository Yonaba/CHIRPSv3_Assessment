root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# Study area map + station distribution along altitude
#   a) DEM (hillshaded) of the study area, rain gauge stations, north arrow,
#      scale bar, 0.05° bounding box and a 3D globe inset (Africa + AOI)
#   b) Number of stations per 50 m altitude interval
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(sf)
  library(terra)
  library(tidyterra)      # geom_spatraster(), scale_fill_hypso_c()
  library(ggnewscale)     # second fill scale (hillshade + DEM)
  library(ggspatial)      # north arrow & scale bar
  library(rnaturalearth)  # world countries (needs rnaturalearthdata)
  library(ggpubr)
  library(fields)         # thin plate spline interpolation (isohyets)
  library(geomtextpath)   # labels drawn along the isohyets
})

# ---- Parameters --------------------------------------------------------------
f_stations  <- "data/csv/p_stations.csv"
f_daily     <- "data/csv/p_stations_daily_data.csv"
f_dem       <- "gis/raster/hydrosheds_dem.tif"
f_aoi       <- "gis/shp/bf_deg.shp"
f_out       <- "figures/0_study_area_stations.png"

period         <- c(1998, 2020)
bbox_buffer    <- 0.05   # degrees, map extent buffer around the study area
display_margin <- 0.10   # degrees, extra room around that extent
alt_bin        <- 50     # m, altitude interval for panel b
target_res     <- 250    # m, DEM resolution used for display (speeds plotting)
z_exag         <- 8      # vertical exaggeration for the hillshade (low relief)
dem_alpha      <- 0.50   # DEM transparency over the hillshade (lower = stronger relief)

inset_frac     <- 0.20   # globe width as a fraction of the map width

# Isohyets
iso_levels     <- seq(300, 1200, 100)   # mm/year
min_frac_year  <- 0.90   # min. fraction of non-missing days for a valid year
min_years      <- 10     # min. number of valid years for a station mean
iso_res        <- 2000   # m, interpolation grid resolution
iso_min_len    <- 30     # km, isohyet pieces shorter than this are dropped

# Legend block (horizontal colour bar + station / isohyet keys), 2°W-2°E
leg_lon   <- c(-2, 2)
leg_top   <- 10.78       # top of the legend background (degrees N)

# Common theme -----------------------------------------------------------------
theme_panel <- theme_bw(base_size = 12) +
  theme(
    plot.title   = element_text(colour = "black", size = 18, hjust = 0, face = "bold"),
    axis.text    = element_text(colour = "black", size = 12),
    axis.title   = element_text(colour = "black", size = 12),
    legend.text  = element_text(colour = "black", size = 12),
    legend.title = element_text(colour = "black", size = 12)
  )

# ---- Data: stations active over 1998-2020 ------------------------------------
daily <- read_csv(f_daily, show_col_types = FALSE) |>
  mutate(Date = mdy(Date)) |>
  filter(between(year(Date), period[1], period[2]))

active <- daily |>
  select(-Date) |>
  summarise(across(everything(), ~ sum(!is.na(.x)))) |>
  pivot_longer(everything(), names_to = "Station", values_to = "n_obs") |>
  mutate(Station = str_trim(Station)) |>
  filter(n_obs > 0)

stations <- read_csv(f_stations, show_col_types = FALSE) |>
  mutate(Station = str_trim(Station)) |>
  filter(Station %in% active$Station) |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

missing <- setdiff(active$Station, stations$Station)
if (length(missing) > 0)
  message("Stations with data but no coordinates: ", paste(missing, collapse = ", "))

# ---- Study area, bounding box and display extent -----------------------------
aoi       <- st_read(f_aoi, quiet = TRUE) |> st_transform(4326) |> st_make_valid()
aoi_union <- st_union(aoi) |> st_make_valid()

bb <- st_bbox(aoi_union)
bbox_ext <- c(xmin = unname(bb["xmin"]) - bbox_buffer, ymin = unname(bb["ymin"]) - bbox_buffer,
              xmax = unname(bb["xmax"]) + bbox_buffer, ymax = unname(bb["ymax"]) + bbox_buffer)
bbox_sf <- st_as_sfc(st_bbox(bbox_ext, crs = st_crs(4326)))

xlim <- unname(c(bbox_ext["xmin"] - display_margin, bbox_ext["xmax"] + display_margin))
ylim <- unname(c(bbox_ext["ymin"] - display_margin, bbox_ext["ymax"] + display_margin))
w <- diff(xlim); h <- diff(ylim)

# ---- DEM + hillshade (computed in EPSG:32630, displayed in EPSG:4326) --------
dem      <- rast(f_dem)
aoi_utm  <- st_transform(aoi_union, 32630)
dem_crop <- crop(dem, vect(st_buffer(aoi_utm, 5000)))

fact  <- max(1, round(target_res / res(dem_crop)[1]))
dem_a <- if (fact > 1) aggregate(dem_crop, fact = fact, fun = "mean", na.rm = TRUE) else dem_crop

slope  <- terrain(dem_a * z_exag, v = "slope",  unit = "radians")
aspect <- terrain(dem_a,          v = "aspect", unit = "radians")
hill   <- shade(slope, aspect, angle = 40, direction = 315, normalize = TRUE)

dem_ll  <- mask(dem_a, vect(aoi_utm)) |> project("EPSG:4326", method = "bilinear")
hill_ll <- mask(hill,  vect(aoi_utm)) |> project("EPSG:4326", method = "bilinear")
names(dem_ll)  <- "elev"
names(hill_ll) <- "hill"
elev_rng <- range(values(dem_ll, mat = FALSE), na.rm = TRUE)

# Station altitudes from the full-resolution DEM
stations$Altitude <- terra::extract(dem, vect(st_transform(stations, 32630)))[, 2]

# ---- Isohyets from station mean annual rainfall (1998-2020) ------------------
annual <- daily |>
  pivot_longer(-Date, names_to = "Station", values_to = "P") |>
  mutate(Station = str_trim(Station),
         P       = suppressWarnings(as.numeric(P)),
         year    = year(Date)) |>
  group_by(Station, year) |>
  summarise(n_ok = sum(!is.na(P)), n_days = n(),
            P = sum(P, na.rm = TRUE), .groups = "drop") |>
  filter(n_ok / n_days >= min_frac_year)

p_mean <- annual |>
  group_by(Station) |>
  summarise(n_years = n(), Pann = mean(P), .groups = "drop") |>
  filter(n_years >= min_years)

st_rain <- stations |>
  inner_join(p_mean, by = "Station") |>
  st_transform(32630)
message("Stations used for isohyets: ", nrow(st_rain),
        " | mean annual rainfall range: ", paste(round(range(st_rain$Pann)), collapse = "-"), " mm")

# Thin plate spline on UTM coordinates, predicted over the study area
tps    <- Tps(st_coordinates(st_rain), st_rain$Pann)
grid_r <- rast(ext(vect(st_buffer(aoi_utm, 10000))), resolution = iso_res, crs = "EPSG:32630")
p_grid <- interpolate(grid_r, tps) |> mask(vect(aoi_utm))

iso <- as.contour(p_grid, levels = iso_levels) |>
  st_as_sf() |>
  st_cast("MULTILINESTRING") |>
  st_cast("LINESTRING") |>
  mutate(len_km = as.numeric(st_length(geometry)) / 1000) |>
  filter(len_km >= iso_min_len) |>
  st_transform(4326) |>
  mutate(label = as.character(level))

# ---- Inset: 3D globe centred on Africa ---------------------------------------
sf_use_s2(TRUE)
lon0 <- 10; lat0 <- 8
R    <- 6371000
ortho_crs <- sprintf("+proj=ortho +lat_0=%s +lon_0=%s +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
                     lat0, lon0)

# Visible hemisphere (spherical cap) used to clip features before projecting
visible <- st_sfc(st_point(c(lon0, lat0)), crs = 4326) |> st_buffer(dist = 9.8e6)

# Natural Earth polygons are valid on the plane but not always on the sphere:
# read them into s2 without checks and repair them with s2_rebuild(). A repaired
# ring can come out with the wrong orientation (the polygon then covers the rest
# of the globe), so polygons larger than a hemisphere are replaced by their
# complement.
world_raw <- ne_countries(scale = 110, returnclass = "sf") |>
  filter(continent != "Antarctica")

world_s2 <- s2::as_s2_geography(st_as_binary(st_geometry(world_raw)), check = FALSE) |>
  s2::s2_rebuild()
flipped <- s2::s2_area(world_s2, radius = R) > 2 * pi * R^2
if (any(flipped))
  world_s2[flipped] <- s2::s2_difference(s2::as_s2_geography(TRUE), world_s2[flipped])
world_s2 <- s2::s2_intersection(world_s2, s2::as_s2_geography(visible))

world <- world_raw |>
  st_set_geometry(st_as_sfc(world_s2)) |>
  filter(!st_is_empty(geometry)) |>
  st_collection_extract("POLYGON") |>
  st_transform(ortho_crs)

africa <- filter(world, continent == "Africa")

grat <- st_graticule(lon = seq(-180, 180, 20), lat = seq(-80, 80, 20)) |>
  st_intersection(visible) |>
  st_transform(ortho_crs)

aoi_glb <- st_transform(aoi_union, ortho_crs)

# Helper: disc in the orthographic plane
disc <- function(x, y, r) st_sfc(st_point(c(x, y)), crs = ortho_crs) |> st_buffer(r, nQuadSegs = 90)

ocean  <- disc(0, 0, R)
halo   <- disc(0, 0, R * 1.05)                 # thin atmosphere
shadow <- disc(R * 0.05, -R * 0.07, R * 1.01)  # drop shadow

# Shading: light source towards the upper left
light <- c(-0.45, 0.45) * R
n_sh  <- 25
t_sh  <- seq_len(n_sh) / n_sh

# Limb darkening: the further from the light, the more layers overlap
dark <- map(t_sh, \(t) st_difference(ocean, disc(light[1] * t, light[2] * t, R * (1 - 0.7 * t)))) |>
  do.call(what = c) |> st_sf(geometry = _)

# Specular highlight around the light point
glow <- map(t_sh, \(t) st_intersection(ocean, disc(light[1], light[2], R * 0.55 * (1 - t) + 1))) |>
  do.call(what = c) |> st_sf(geometry = _)

p_globe <- ggplot() +
  geom_sf(data = shadow,  fill = "grey30", alpha = 0.25, colour = NA) +
  geom_sf(data = halo,    fill = "#cfe6fb", alpha = 0.8, colour = NA) +
  geom_sf(data = ocean,   fill = "#8fbfe0", colour = NA) +
  geom_sf(data = grat,    colour = "white", alpha = 0.5, linewidth = 0.15) +
  geom_sf(data = world,   fill = "grey88", colour = "grey60", linewidth = 0.1) +
  geom_sf(data = africa,  fill = "#cdbf9f", colour = "grey45", linewidth = 0.1) +
  geom_sf(data = aoi_glb, fill = "darkorange2", colour = "darkorange4", linewidth = 0.25) +
  geom_sf(data = dark,    fill = "black", alpha = 0.035, colour = NA) +
  geom_sf(data = glow,    fill = "white", alpha = 0.03,  colour = NA) +
  geom_sf(data = ocean,   fill = NA, colour = "grey25", linewidth = 0.3) +
  coord_sf(crs = ortho_crs, expand = FALSE,
           xlim = c(-1.06, 1.08) * R, ylim = c(-1.10, 1.06) * R) +
  theme_void()

# Inset position (top-left corner, in map coordinates)
inset_size <- inset_frac * w
inset_xmin <- xlim[1] + 0.005 * w
inset_ymax <- ylim[2] - 0.005 * h

# ---- Legend block drawn in map coordinates -----------------------------------
# Drawn manually so it spans exactly leg_lon whatever the export size.
# Vertical layout (degrees below leg_top):
leg_box  <- c(leg_top - 1.06, leg_top)            # background bottom / top
title_y  <- leg_top - 0.20                        # "Elevation (m)" baseline
bar_y    <- leg_top - c(0.46, 0.30)               # colour bar bottom / top
key_y    <- leg_top - 0.88                        # station / isohyet keys row

n_col    <- 200
bar_x    <- seq(leg_lon[1], leg_lon[2], length.out = n_col + 1)
bar_cols <- hypso.colors2(n_col, palette = "dem_poster")

brk   <- pretty(elev_rng)
brk   <- brk[brk >= elev_rng[1] & brk <= elev_rng[2]]
brk_x <- leg_lon[1] + diff(leg_lon) * (brk - elev_rng[1]) / diff(elev_rng)

txt_size <- 12 / .pt   # matches size 12 text elements
iso_col  <- "blue"

# ---- Panel a: main map -------------------------------------------------------
p_map <- ggplot() +
  geom_spatraster(data = hill_ll, show.legend = FALSE) +
  scale_fill_gradient(low = "grey5", high = "white", na.value = "transparent", guide = "none") +
  new_scale_fill() +
  geom_spatraster(data = dem_ll, alpha = dem_alpha) +
  scale_fill_hypso_c(palette = "dem_poster", limits = elev_rng,
                     na.value = "transparent", guide = "none") +
  geom_sf(data = aoi,       fill = NA, colour = "grey35", linewidth = 0.25) +
  geom_sf(data = aoi_union, fill = NA, colour = "black",  linewidth = 0.6) +
  # Isohyets, labelled along the line (mm)
  geom_textsf(data = iso, aes(label = label), colour = iso_col, linecolour = iso_col,
              linewidth = 0.6, size = 3.5, fontface = "bold", hjust = 0.5,
              remove_long = TRUE) +
  geom_sf(data = stations,  colour = "black", size = 1.8, shape = 16) +
  # Legend background
  annotate("rect", xmin = leg_lon[1] - 0.2, xmax = leg_lon[2] + 0.2,
           ymin = leg_box[1], ymax = leg_box[2],
           fill = alpha("white", 0.75), colour = NA) +
  # Elevation colour bar
  annotate("text", x = leg_lon[1], y = title_y, label = "Elevation (m)",
           hjust = 0, vjust = 0, size = txt_size, colour = "black") +
  annotate("rect", xmin = head(bar_x, -1), xmax = tail(bar_x, -1),
           ymin = bar_y[1], ymax = bar_y[2], fill = bar_cols, colour = NA) +
  annotate("rect", xmin = leg_lon[1], xmax = leg_lon[2],
           ymin = bar_y[1], ymax = bar_y[2], fill = NA, colour = "black", linewidth = 0.3) +
  annotate("segment", x = brk_x, xend = brk_x,
           y = bar_y[1], yend = bar_y[1] - 0.04, colour = "black", linewidth = 0.3) +
  annotate("text", x = brk_x, y = bar_y[1] - 0.07, label = brk,
           vjust = 1, size = txt_size, colour = "black") +
  # Station key
  annotate("point", x = leg_lon[1] + 0.08, y = key_y, colour = "black", size = 2.2, shape = 16) +
  annotate("text", x = leg_lon[1] + 0.20, y = key_y,
           label = sprintf("Rainfall station (%d)", nrow(stations)),
           hjust = 0, vjust = 0.5, size = txt_size, colour = "black") +
  # Isohyet key
  annotate("segment", x = 0.30, xend = 0.75, y = key_y, yend = key_y,
           colour = iso_col, linewidth = 0.8) +
  annotate("text", x = 0.85, y = key_y, label = "Annual rainfall (mm)",
           hjust = 0, vjust = 0.5, size = txt_size, colour = "black") +
  annotation_custom(ggplotGrob(p_globe),
                    xmin = inset_xmin, xmax = inset_xmin + inset_size,
                    ymin = inset_ymax - inset_size, ymax = inset_ymax) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1.5, "cm"), width = unit(1.5, "cm"),
                         pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"),
                         style = north_arrow_fancy_orienteering(text_col = "black")) +
  annotation_scale(location = "br", unit_category = "metric", width_hint = 0.25,
                   text_cex = 1, pad_x = unit(0.4, "cm"), pad_y = unit(0.4, "cm")) +
  coord_sf(crs = 4326, xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "a) Study area and rain gauge stations", x = NULL, y = NULL) +
  theme_panel

# ---- Panel b: station distribution along altitude ----------------------------
alt_df <- stations |>
  st_drop_geometry() |>
  filter(!is.na(Altitude)) |>
  mutate(bin_low = floor(Altitude / alt_bin) * alt_bin) |>
  count(bin_low, name = "n")

alt_df <- alt_df |>
  complete(bin_low = seq(min(bin_low) - alt_bin, max(bin_low) + alt_bin, by = alt_bin),
           fill = list(n = 0)) |>
  mutate(alt_mid = bin_low + alt_bin / 2) |>
  arrange(alt_mid)

p_alt <- ggplot(alt_df, aes(x = n, y = alt_mid)) +
  geom_path(colour = "blue", linewidth = 1, linetype = "solid") +
  scale_y_continuous(breaks = seq(min(alt_df$bin_low), max(alt_df$bin_low) + alt_bin, by = alt_bin)) +
  scale_x_continuous(breaks = function(l) unique(floor(pretty(l))),
                     limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "b) Stations by altitude",
       x = "Number of stations", y = "Altitude (m)") +
  theme_panel

# ---- Assemble & export -------------------------------------------------------
right_col <- ggarrange(p_alt, ncol = 1)   # add further secondary panels here
fig <- ggarrange(p_map, right_col, ncol = 2, widths = c(2.3, 1))

dir.create("figures", showWarnings = FALSE)
ggsave(f_out, fig, width = 16, height = 8.5, units = "in",
       dpi = 500, scale = 0.8, bg = "white")