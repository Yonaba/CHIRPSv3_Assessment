root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# -----------------------------------------------------------------------------
# Number of stations used by CHIRPS v2 and CHIRPS v3 over Burkina Faso,
# per 1° latitude band (9-16°N), monthly, 1998-01 to 2020-12.
# Output: figures/0_cmp_chirps_stations_used.png
# -----------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(ggpubr)

# ---- Settings ---------------------------------------------------------------
files <- c(v2 = file.path("data", "chirps_stations", "chirps_stations_v2.csv"),
           v3 = file.path("data", "chirps_stations", "chirps_stations_v3.csv"))
out_png   <- file.path("figures", "0_cmp_chirps_stations_used.png")
bands      <- c(9, 11, 13)               # lower edge of each band (9-11, 11-13, 13-15)
band_width <- 2                          # band width (°)
ylim_band  <- c(0, 90)                   # y-axis limits, latitude band panels
ylim_total <- c(0, 150)                  # y-axis limits, all-bands panel
smooth_span <- 0.1                       # loess span for the smoothed line
date_from <- as.Date("1998-01-01")
date_to   <- as.Date("2020-12-01")

dir.create("figures", showWarnings = FALSE)

# ---- Load data --------------------------------------------------------------
load_stations <- function(path) {
  dt <- fread(path, showProgress = FALSE)
  
  # Locate the latitude / longitude columns whatever their exact name
  lat_col <- grep("^lat", names(dt), ignore.case = TRUE, value = TRUE)[1]
  lon_col <- grep("^lon", names(dt), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(lat_col)) stop("No latitude column found in ", path)
  
  dt <- dt[, .(year, month,
               lat = as.numeric(get(lat_col)),
               lon = if (!is.na(lon_col)) as.numeric(get(lon_col)) else NA_real_)]
  
  # One station = one location per month (guards against duplicates,
  # e.g. a station listed in both the v2 "global" and "extra" files)
  unique(dt, by = c("year", "month", "lat", "lon"))
}

stations <- rbindlist(lapply(files, load_stations), idcol = "version")

# ---- Count stations per latitude band and month -----------------------------
stations <- stations[lat >= min(bands) & lat < max(bands) + band_width]
stations[, band := min(bands) + band_width * floor((lat - min(bands)) / band_width)]
stations[, date := as.Date(sprintf("%d-%02d-01", year, month))]
stations <- stations[date >= date_from & date <= date_to]

counts <- stations[, .(n = .N), by = .(version, band, date)]

# Complete grid so months without any station show as 0, not as a gap
full <- CJ(version = names(files), band = bands,
           date = seq(date_from, date_to, by = "month"))
counts <- counts[full, on = .(version, band, date)]
counts[is.na(n), n := 0L]

totals <- counts[, .(n = sum(n)), by = .(version, date)]

# ---- Plot helpers -----------------------------------------------------------
base_theme <- theme_bw() +
  theme(
    axis.title  = element_text(size = 12, colour = "black"),
    axis.text   = element_text(size = 12, colour = "black"),
    plot.title  = element_text(size = 18, face = "bold", hjust = 0, colour = "black"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(3, 8, 3, 5)
  )

x_scale <- scale_x_date(limits = c(date_from, date_to),
                        date_breaks = "2 years", date_labels = "%Y",
                        expand = c(0, 0))

make_panel <- function(dt, label, ylim, left_col, bottom) {
  p <- ggplot(dt, aes(x = date, y = n)) +
    geom_line(colour = "black", linewidth = 0.3, linetype = "dotted") +
    geom_smooth(method = "loess", formula = y ~ x, span = smooth_span,
                se = FALSE, colour = "darkorange3", linewidth = 0.9) +
    annotate("text", x = date_from, y = Inf, label = label,
             hjust = -0.05, vjust = 1.4, size = 12 / .pt,
             colour = "black") +
    x_scale +
    # coord_cartesian keeps the loess fit computed on all data, even where it
    # would dip slightly below 0
    coord_cartesian(ylim = ylim, expand = FALSE) +
    labs(x = if (bottom) "" else NULL,
         y = if (left_col) "No. of stations" else NULL) +
    base_theme
  
  if (!bottom)   p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  if (!left_col) p <- p + theme(axis.title.y = element_blank())
  p
}

make_column <- function(ver, title, left_col) {
  # Band panels, northernmost (13-15°N) on top
  band_panels <- lapply(rev(bands), function(b) {
    make_panel(counts[version == ver & band == b],
               label    = sprintf("%d\u2013%d\u00b0N", b, b + band_width),
               ylim     = ylim_band,
               left_col = left_col,
               bottom   = FALSE)
  })
  
  total_panel <- make_panel(totals[version == ver],
                            label    = sprintf("All bands (%d\u2013%d\u00b0N)",
                                               min(bands), max(bands) + band_width),
                            ylim     = ylim_total,
                            left_col = left_col,
                            bottom   = TRUE)
  
  col <- ggarrange(plotlist = c(band_panels, list(total_panel)),
                   ncol = 1, align = "v",
                   heights = c(rep(1, length(bands)), 2))
  
  annotate_figure(col, top = text_grob(title, size = 16, face = "bold",
                                       color = "black", hjust = 0, x = 0.01))
}

# ---- Assemble and export ----------------------------------------------------
col_v2 <- make_column("v2", "a) Stations in CHIRPS v2\n", left_col = TRUE)
col_v3 <- make_column("v3", "b) Stations in CHIRPS v3\n", left_col = FALSE)

fig <- ggarrange(col_v2, col_v3, ncol = 2, widths = c(1.05, 1))

ggsave(out_png, fig, width = 18, height = 11, units = "in",
       dpi = 350, scale = 0.8, bg = "white")

message("Figure saved to ", out_png)
