root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# Maps of station-level KGE (CHIRPS v2) and KGE differences (v3 - v2) for CHIRPS v3 (RNL)
# and CHIRPS v3 (SAT) against station observations, 1998-2020, at 5 timescales
# (daily, pentad, dekadal, monthly, annual).
# Output: multi-panel figure (5 rows = timescales x 3 columns = products)
#         in "figures/2_chirps_kge_diff_maps.png"
# =============================================================================

# ---- 0. Packages -------------------------------------------------------------
pkgs <- c("terra", "sf", "dplyr", "tidyr", "lubridate", "ggplot2",
          "ggspatial", "ggpubr", "readr")
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
outline_file  <- "gis/shp/bf_deg.shp"

products <- c("CHIRPS v2"       = "chirps_v2",
              "CHIRPS v3 (RNL)" = "chirps_v3_rnl",
              "CHIRPS v3 (SAT)" = "chirps_v3_sat")

timescales <- c("Daily", "Pentad", "Dekadal", "Monthly", "Annual")

# Minimum number of paired (complete) values required to compute KGE of a
# station at a given timescale (about one year of data; 72 pentads and
# 36 dekads per calendar year)
min_pairs <- c(Daily = 365, Pentad = 72, Dekadal = 36, Monthly = 24, Annual = 8)

# Colour scale limits for KGE (values outside are squished to the end colours)
kge_limits <- c(-1, 1)

dir.create("figures", showWarnings = FALSE)

# ---- 2. Load spatial data and station locations -----------------------------
outline <- st_read(outline_file, quiet = TRUE) |> st_transform(4326)

bb <- st_bbox(outline)
xlim_map <- c(bb[["xmin"]] - 0.05, bb[["xmax"]] + 0.05)
ylim_map <- c(bb[["ymin"]] - 0.05, bb[["ymax"]] + 0.05)

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

common_st <- intersect(stations$Station, unique(obs_long$Station))
if (length(common_st) == 0) stop("No station names match between the two csv files.")
if (length(common_st) < nrow(stations))
  message(nrow(stations) - length(common_st),
          " station(s) from p_stations.csv have no data column and are dropped.")
stations <- stations |> filter(Station %in% common_st)

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

# ---- 5. KGE (Kling et al. 2012) ---------------------------------------------
calc_kge <- function(obs, sim) {
  ok <- is.finite(obs) & is.finite(sim)
  o <- obs[ok]; s <- sim[ok]
  if (length(o) < 3) return(NA_real_)
  if (!(sd(o) > 0 && sd(s) > 0 && mean(o) > 0 && mean(s) > 0)) return(NA_real_)
  r     <- cor(o, s, method = "pearson")
  beta  <- mean(s) / mean(o)
  gamma <- (sd(s) / mean(s)) / (sd(o) / mean(o))
  1 - sqrt((r - 1)^2 + (beta - 1)^2 + (gamma - 1)^2)
}

# Accumulation key for each timescale (CHIRPS calendar-month based periods)
# - Pentad : 6 per month (1-5, 6-10, 11-15, 16-20, 21-25, 26-end of month)
# - Dekadal: 3 per month (1-10, 11-20, 21-end of month)
add_period <- function(df, scale) {
  switch(scale,
         Daily   = df |> mutate(period = as.character(Date)),
         Pentad  = df |> mutate(period = paste(year(Date), month(Date),
                                               pmin(ceiling(mday(Date) / 5), 6))),
         Dekadal = df |> mutate(period = paste(year(Date), month(Date),
                                               pmin(ceiling(mday(Date) / 10), 3))),
         Monthly = df |> mutate(period = paste(year(Date), month(Date))),
         Annual  = df |> mutate(period = as.character(year(Date))))
}

# A period is valid only if all its days are valid (paired obs & sim)
accumulate <- function(df, scale) {
  if (scale == "Daily") return(df |> select(Station, obs, sim))
  df |>
    add_period(scale) |>
    group_by(Station, period) |>
    summarise(obs = if (anyNA(obs)) NA_real_ else sum(obs),
              sim = if (anyNA(sim)) NA_real_ else sum(sim),
              .groups = "drop")
}

kge_all <- list()
for (p in names(products)) {
  message("Computing KGE for ", p, " ...")
  paired <- obs_long |>
    inner_join(sim_list[[p]], by = c("Date", "Station")) |>
    mutate(bad = is.na(obs) | is.na(sim),
           obs = ifelse(bad, NA_real_, obs),
           sim = ifelse(bad, NA_real_, sim)) |>
    select(-bad)
  
  for (sc in timescales) {
    agg <- accumulate(paired, sc)
    kge_all[[paste(p, sc)]] <- agg |>
      group_by(Station) |>
      summarise(n   = sum(is.finite(obs) & is.finite(sim)),
                KGE = calc_kge(obs, sim),
                .groups = "drop") |>
      mutate(KGE = ifelse(n >= min_pairs[[sc]], KGE, NA_real_),
             Source = p, Timescale = sc)
  }
}

kge <- bind_rows(kge_all) |>
  left_join(stations |> select(Station, Longitude, Latitude), by = "Station") |>
  mutate(Source    = factor(Source, levels = names(products)),
         Timescale = factor(Timescale, levels = timescales))

write_csv(kge, "data/csv/2_station_kge_long.csv")

# ---- 6. Differences with respect to CHIRPS v2 -------------------------------
# Column 1: KGE of CHIRPS v2. Columns 2-3: KGE(v3) - KGE(v2) per station
# (only stations with a valid KGE for both products).
ref <- "CHIRPS v2"

kge_wide <- kge |>
  select(Station, Longitude, Latitude, Timescale, Source, KGE) |>
  pivot_wider(names_from = Source, values_from = KGE)

delta_tbl <- bind_rows(lapply(setdiff(names(products), ref), function(pr) {
  kge_wide |>
    transmute(Station, Longitude, Latitude, Timescale,
              Source = pr, delta_KGE = .data[[pr]] - .data[[ref]])
}))
write_csv(delta_tbl, "data/csv/2_station_delta_kge_long.csv")

panel_data <- function(sc, pr) {
  base <- kge_wide |> filter(Timescale == sc)
  if (pr == ref) {
    base |> transmute(Station, Longitude, Latitude, value = .data[[ref]])
  } else {
    base |> transmute(Station, Longitude, Latitude,
                      value = .data[[pr]] - .data[[ref]])
  }
}

# ---- 7. Maps -----------------------------------------------------------------
delta_limits <- c(-0.3, 0.3)   # symmetric colour range for the differences
delta_tol    <- 0.02           # |change| below this is drawn with a smaller dot

n_row <- length(timescales)
n_col <- length(products)

plots <- list()
k <- 0
for (i in seq_along(timescales)) {          # rows
  for (j in seq_along(products)) {          # columns
    k <- k + 1
    sc <- timescales[i]
    pr <- names(products)[j]
    is_delta <- pr != ref
    
    d <- panel_data(sc, pr) |>
      filter(is.finite(value)) |>
      mutate(small = if (is_delta) abs(value) < delta_tol else FALSE) |>
      arrange(abs(value))                   # largest values plotted on top
    
    ttl <- if (is_delta) {
      paste0(letters[k], ") ΔKGE (", sub("CHIRPS ", "", pr), " − v2) - ", sc)
    } else {
      paste0(letters[k], ") ", pr, " - ", sc)
    }
    
    fill_scale <- if (is_delta) {
      scale_fill_gradientn(
        colours = c("red3", "white", "green4"),
        limits  = delta_limits, oob = scales::squish,
        name    = "ΔKGE (v3 − v2)",
        guide   = guide_colourbar(direction = "horizontal",
                                  title.position = "left", title.vjust = 0.8,
                                  barwidth = unit(10, "cm"),
                                  barheight = unit(0.6, "cm")))
    } else {
      scale_fill_gradientn(
        colours = c("red3", "yellow", "green4"),
        limits  = kge_limits, oob = scales::squish,
        name    = "KGE (CHIRPS v2)",
        guide   = guide_colourbar(direction = "horizontal",
                                  title.position = "left", title.vjust = 0.8,
                                  barwidth = unit(6, "cm"),
                                  barheight = unit(0.6, "cm")))
    }
    
    p <- ggplot() +
      geom_sf(data = outline, fill = "grey97", colour = "black", linewidth = 0.4) +
      geom_point(data = d, aes(x = Longitude, y = Latitude, fill = value, size = small),
                 shape = 21, colour = "black", stroke = 0.3) +
      scale_size_manual(values = c(`FALSE` = 3.2, `TRUE` = 2), guide = "none") +
      fill_scale +
      coord_sf(xlim = xlim_map, ylim = ylim_map, expand = FALSE) +
      annotation_north_arrow(location = "tr", which_north = "true",
                             style = north_arrow_minimal(text_size = 10),
                             height = unit(1.2, "cm"), width = unit(1.2, "cm")) +
      annotation_scale(location = "br", unit_category = "metric",
                       width_hint = 0.3, text_cex = 0.9) +
      labs(title = ttl, x = "Longitude", y = "Latitude") +
      theme_bw() +
      theme(plot.title   = element_text(colour = "black", size = 18,
                                        face = "bold", hjust = 0),
            axis.text    = element_text(colour = "black", size = 12),
            axis.title   = element_text(colour = "black", size = 12),
            legend.text  = element_text(colour = "black", size = 12),
            legend.title = element_text(colour = "black", size = 12),
            legend.position = "bottom",
            panel.grid.minor = element_blank())
    
    # x-axis title only on bottom row, y-axis title only on leftmost column
    if (i != n_row) p <- p + theme(axis.title.x = element_blank())
    if (j != 1)     p <- p + theme(axis.title.y = element_blank())
    
    plots[[k]] <- p
  }
}

# ---- 8. Assemble: 2 legends (KGE for v2, delta KGE for v3 maps) --------------
leg_kge   <- ggpubr::as_ggplot(ggpubr::get_legend(plots[[1]]))   # v2 KGE
leg_delta <- ggpubr::as_ggplot(ggpubr::get_legend(plots[[2]]))   # delta KGE

plots_nl <- lapply(plots, function(p) p + theme(legend.position = "none"))

grid_maps <- ggarrange(plotlist = plots_nl, ncol = n_col, nrow = n_row, align = "hv")
legends   <- ggarrange(leg_kge, leg_delta, ncol = 2, widths = c(1, 2))
fig <- ggarrange(grid_maps, legends, ncol = 1, heights = c(1, 0.04))

ggsave("figures/2_chirps_kge_diff_maps.png", fig,
       width = 22, height = 26, units = "in", dpi = 500, bg = "white",
       scale = 0.8)