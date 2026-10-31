root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# ------------------------------------------------------------------------------
# Half Taylor diagrams: CHIRPS products (v2, v3 RNL, v3 SAT) vs observed daily
# precipitation (pooled over all stations), at 5 timescales:
# daily, pentad, dekadal, monthly, annual.
# One diagram per timescale, arranged in a single figure with ggpubr::ggarrange.
# A table of performance metrics (MAE, RMSE, PBIAS, R2, NSE, KGE-2012) per
# product and timescale is also exported.
#
# Aggregation rules (done station by station, BEFORE pooling):
#  - pentad: CHIRPS calendar pentads, 6 per month (1-5, 6-10, 11-15, 16-20,
#    21-25, 26-end of month) -> 72 pentads per year. The 6th pentad has 3 to
#    6 days depending on the month.
#  - dekadal: CHIRPS calendar dekads, 3 per month (1-10, 11-20, 21-end of
#    month) -> 36 dekads per year. The 3rd dekad has 8 to 11 days.
#  - monthly / annual: calendar months / calendar years.
#  - A sum is NA as soon as one day of the block is missing.
# Only (obs, all products) pairs that are complete are kept, so that all
# products are evaluated on exactly the same observations.
# ------------------------------------------------------------------------------

# remotes::install_github("AlexandreWadoux/modelskill")
library(modelskill)
library(ggplot2)
library(ggpubr)
library(hydroGOF)

# ---- Settings ----------------------------------------------------------------
csv_folder   <- file.path("data", "csv")
fig_folder   <- "figures"
table_folder <- "tables"

obs_file <- file.path(csv_folder, "p_stations_daily_data.csv")

# c("name shown in the legend" = "file")
products <- c(
  "CHIRPS v2"       = file.path(csv_folder, "p_chirpsv2_daily_data.csv"),
  "CHIRPS v3 (RNL)" = file.path(csv_folder, "p_chirpsv3rnl_daily_data.csv"),
  "CHIRPS v3 (SAT)" = file.path(csv_folder, "p_chirpsv3sat_daily_data.csv")
)

# Legend order and dot colours (names must match the legend labels)
obs_label   <- "Observations"
dot_colours <- c(
  "Observations"    = "black",
  "CHIRPS v2"       = "darkred",
  "CHIRPS v3 (RNL)" = "darkblue",
  "CHIRPS v3 (SAT)" = "darkgreen"
)
dot_size <- 4

# Timescales (order of the panels: a), b), c), ...)
timescales <- c("daily", "pentad", "dekadal", "monthly", "annual")
ts_labels  <- c(
  daily   = "Daily",
  pentad  = "Pentad (5 days)",
  dekadal = "Dekadal (10 days)",
  monthly = "Monthly",
  annual  = "Annual"
)

# Common period
start_date <- as.Date("1998-01-01")
end_date   <- as.Date("2020-12-31")

dir.create(fig_folder,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_folder, showWarnings = FALSE, recursive = TRUE)

# ---- Helpers -----------------------------------------------------------------
# Wide daily file: first column = Date ("m/d/Y"), then one column per station
read_wide <- function(file) {
  df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  names(df)[1]  <- "Date"
  names(df)[-1] <- trimws(names(df)[-1])
  df$Date <- as.Date(trimws(as.character(df$Date)), format = "%m/%d/%Y")
  if (anyNA(df$Date)) stop("Some dates could not be parsed (expected m/d/Y) in ", file)
  df
}

# ---- Load data ---------------------------------------------------------------
message("Reading data ...")
all_dates <- seq(start_date, end_date, by = "day")
raw <- c(list(obs = read_wide(obs_file)), lapply(products, read_wide))
names(raw) <- c(obs_label, names(products))

# Common stations across observations and all products
stations <- Reduce(intersect, lapply(raw, function(x) setdiff(names(x), "Date")))
if (length(stations) == 0) stop("No common station between obs and products.")
n_max <- max(sapply(raw, function(x) length(setdiff(names(x), "Date"))))
if (length(stations) < n_max) {
  warning("Only ", length(stations), " stations are common to obs and all products.")
}

# Daily matrices (days x stations) on the complete 1998-2020 date sequence
daily <- lapply(raw, function(x) {
  m <- as.matrix(x[match(all_dates, x$Date), stations, drop = FALSE])
  storage.mode(m) <- "double"
  m
})

# ---- Aggregation to the different timescales ----------------------------------
yr  <- as.integer(format(all_dates, "%Y"))
mon <- as.integer(format(all_dates, "%m"))
dom <- as.integer(format(all_dates, "%d"))

# Calendar-month based CHIRPS periods
pentad_in_month <- pmin((dom - 1) %/% 5  + 1, 6)   # 1..6  (26-end -> 6)
dekad_in_month  <- pmin((dom - 1) %/% 10 + 1, 3)   # 1..3  (21-end -> 3)

# For each timescale: grouping key per day + number of days required in a group
# (NA = keep all groups). Since all_dates covers complete calendar years, every
# group is complete by construction; missing days are handled by NA propagation.
groups <- list(
  daily   = list(key = as.character(all_dates),            need = 1),
  pentad  = list(key = paste(yr, mon, pentad_in_month),    need = NA),
  dekadal = list(key = paste(yr, mon, dekad_in_month),     need = NA),
  monthly = list(key = format(all_dates, "%Y-%m"),         need = NA),
  annual  = list(key = as.character(yr),                   need = NA)
)

# Sum per group and per station. NA is propagated (rowsum with na.rm = FALSE).
aggregate_sum <- function(m, key, need) {
  f    <- factor(key, levels = unique(key))
  s    <- rowsum(m, f, reorder = FALSE)
  cnt  <- as.integer(table(f))
  keep <- if (is.na(need)) rep(TRUE, length(cnt)) else cnt == need
  s[keep, , drop = FALSE]
}

# Pooled vectors (all stations stacked) for one timescale
pool_timescale <- function(ts) {
  g   <- groups[[ts]]
  agg <- lapply(daily, aggregate_sum, key = g$key, need = g$need)
  pooled <- lapply(agg, function(m) as.numeric(m))   # stack stations column by column
  
  # Keep only complete obs/sim pairs (same pairs for every product)
  ok <- complete.cases(do.call(cbind, pooled))
  pooled <- lapply(pooled, function(x) x[ok])
  list(obs = pooled[[obs_label]], sims = pooled[names(products)])
}

pooled_all <- lapply(timescales, function(ts) {
  message("Aggregating: ", ts)
  pool_timescale(ts)
})
names(pooled_all) <- timescales

# ---- Performance table ---------------------------------------------------------
# hydroGOF convention: PBIAS = 100 * sum(sim - obs) / sum(obs)
# PBIAS > 0 = overestimation, < 0 = underestimation
metrics_row <- function(sim, obs) {
  data.frame(
    n    = length(obs),
    MAE  = hydroGOF::mae(sim, obs),
    RMSE = hydroGOF::rmse(sim, obs),
    PBIAS = hydroGOF::pbias(sim, obs),
    R2   = hydroGOF::R2(sim, obs),
    NSE  = hydroGOF::NSE(sim, obs),
    KGE  = hydroGOF::KGE(sim, obs, method = "2012")
  )
}

tab <- do.call(rbind, lapply(timescales, function(ts) {
  p <- pooled_all[[ts]]
  do.call(rbind, lapply(names(products), function(prod) {
    cbind(Timescale = ts_labels[[ts]], Product = prod,
          metrics_row(p$sims[[prod]], p$obs))
  }))
}))
rownames(tab) <- NULL
num_cols <- setdiff(names(tab), c("Timescale", "Product", "n"))
tab[num_cols] <- lapply(tab[num_cols], round, 4)

write.csv(tab, file.path(table_folder, "1_performance_metrics_precip.csv"),
          row.names = FALSE)

# ---- Theme (same layout as the template) ---------------------------------------
theme_panel <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title   = element_text(hjust = 0, face = "bold", colour = "black", size = 18),
      axis.title   = element_text(colour = "black", size = 12),
      axis.text    = element_text(colour = "black", size = 12),
      legend.position      = "bottom",
      legend.direction     = "horizontal",
      legend.justification = "center",
      legend.title = element_blank(),
      legend.text  = element_text(colour = "black", size = 12)
    )
}

# ---- Build one Taylor diagram per timescale ------------------------------------
fig_ncol <- 3
panel_letters <- paste0(letters[seq_along(timescales)], ")")

# Observations sit at the Taylor reference point (sigma* = 1, r = 1)
obs_point <- data.frame(x = 1, y = 0, Model = obs_label)

plots <- lapply(seq_along(timescales), function(i) {
  ts <- timescales[i]
  message("Plotting: ", ts)
  pool <- pooled_all[[ts]]
  
  p <- gg_taylor(
    mods       = pool$sims,
    obs        = pool$obs,
    half       = TRUE,
    legend     = TRUE,
    point_size = dot_size
  ) +
    geom_point(data = obs_point, aes(x = x, y = y, colour = Model),
               size = dot_size) +
    scale_colour_manual(values = dot_colours, breaks = names(dot_colours),
                        name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = dot_size))) +
    labs(title = paste(panel_letters[i], ts_labels[[ts]])) +
    ylab("") +
    theme_panel()
  
  # x-axis title transparent (not removed) in the first row so that all panels
  # keep the same dimensions
  if (i <= fig_ncol) p <- p + theme(axis.title.x = element_text(colour = "transparent"))
  p
})
names(plots) <- timescales

# ---- Arrange and export --------------------------------------------------------
fig <- ggarrange(
  plotlist      = plots,
  ncol          = fig_ncol,
  nrow          = 2,
  common.legend = TRUE,
  legend        = "bottom",
  align         = "v"
)

fig <- fig + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(
  filename = file.path(fig_folder, "1_taylor_precip_timescales.png"),
  plot     = fig,
  width    = 13,
  height   = 9,
  units    = "in",
  dpi      = 350,
  bg       = "white"
)

message("Done.")