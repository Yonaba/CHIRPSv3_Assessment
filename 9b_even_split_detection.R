root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# 9b - Detection of the even-split rule in the daily disaggregation
#      Observations vs CHIRPS v2, CHIRPS v3 (RNL) and CHIRPS v3 (SAT), 1998-2020
# -----------------------------------------------------------------------------
# Rationale
#   When the daily driver has no rain in a pentad that CHIRPS considers wet, the
#   pentad total can only be spread evenly. Such even splits create artificial
#   light rain (and artificial early-season wet days). Gauges almost never
#   produce exactly equal daily amounts, so they give the "chance" baseline.
#
# Method
#   1. Same daily sample and complete CHIRPS pentads as in 9a (6 pentads per
#      month, 3 to 6 days). Values < PMIN are set to 0.
#   2. Each station-pentad with a non-zero total is classified as:
#        "Exact 1/n split"        all n days wet and equal
#        "Equal split, wet days"  2 to n-1 wet days, all equal
#        "Single wet day"         one wet day only (not an even split)
#        "Uneven"                 anything else
#      "Equal" = max - min of the non-zero values <= max(tol_abs, tol_rel * mean).
#   3. Even-split frequency (% of non-zero pentads) by series, month, zone and
#      pentad-total class; share of the wet days of each intensity class that
#      come from even-split pentads.
#
# Output: figures/9b_even_split_detection.png
# =============================================================================

# ---- 0. Packages -------------------------------------------------------------
pkgs <- c("data.table", "ggplot2", "patchwork")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 1. Settings -------------------------------------------------------------
start_date <- as.Date("1998-01-01")
end_date   <- as.Date("2020-12-31")
all_dates  <- seq(start_date, end_date, by = "day")

f_stations <- "data/csv/p_stations.csv"
f_obs      <- "data/csv/p_stations_daily_data.csv"
f_prod     <- c(V2  = "data/csv/p_chirpsv2_daily_data.csv",
                RNL = "data/csv/p_chirpsv3rnl_daily_data.csv",
                SAT = "data/csv/p_chirpsv3sat_daily_data.csv")
fig_file   <- "figures/9b_even_split_detection.png"

PMIN         <- 0.1             # mm: drizzle / wet-day threshold
tol_abs      <- 0.01            # mm: absolute tolerance for "equal" daily values
tol_rel      <- 0.001           # relative tolerance (fraction of the mean wet-day value)
zone_method  <- "isohyet"       # "isohyet", "latitude" or "column" (see 9a)
isohyets     <- c(600, 900)
lat_breaks   <- c(11.5, 14)
min_cov_year <- 0.90
min_n        <- 30              # min. number of non-zero pentads to plot a month/zone value
plot_months  <- 6:10            # months shown in panel a

tot_breaks <- c(0, 5, 15, 30, 60, Inf)
tot_labels <- c("0-5", "5-15", "15-30", "30-60", ">60")
int_breaks <- c(PMIN, 1, 2, 5, 10, 20, Inf)       # daily intensity classes (mm)
int_labels <- c("0.1-1", "1-2", "2-5", "5-10", "10-20", ">20")

src_lab  <- c(OBS = "Observations", V2 = "CHIRPS v2",
              RNL = "CHIRPS v3 (RNL)", SAT = "CHIRPS v3 (SAT)")
src_cols <- c("Observations"    = "grey45",
              "CHIRPS v2"       = "darkred",
              "CHIRPS v3 (RNL)" = "darkblue",
              "CHIRPS v3 (SAT)" = "darkgreen")
type_lab  <- c(full = "Exact 1/n split (all days)", partial = "Equal split across wet days")
type_cols <- setNames(c("#8c2d04", "#fdae6b"), type_lab)
zone_levels <- c("Sahelian", "Sudano-Sahelian", "Sudanian")

dir.create("figures", showWarnings = FALSE)

# ---- 2. Data -----------------------------------------------------------------
read_daily <- function(f) {
  if (!file.exists(f)) stop(f, " not found (run 1_extract_stations_daily.R first).")
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  names(d) <- trimws(names(d)); names(d)[1] <- "Date"
  dts <- as.Date(as.character(d$Date), format = "%m/%d/%Y")
  x <- d[match(all_dates, dts), -1, drop = FALSE]
  m <- matrix(suppressWarnings(as.numeric(unlist(x, use.names = FALSE))),
              nrow = nrow(x), dimnames = list(NULL, names(x)))
  m[!is.na(m) & m < 0]    <- NA
  m[!is.na(m) & m < PMIN] <- 0
  m
}

stn <- read.csv(f_stations, check.names = FALSE, stringsAsFactors = FALSE)
stn$Station <- trimws(as.character(stn$Station))
stn <- stn[!is.na(stn$Longitude) & !is.na(stn$Latitude) & !duplicated(stn$Station), ]

obs_m    <- read_daily(f_obs)
st_names <- intersect(stn$Station, colnames(obs_m))
stn      <- stn[match(st_names, stn$Station), ]
n_st     <- length(st_names)
message(n_st, " stations used.")

mats <- list(OBS = obs_m[, st_names, drop = FALSE])
for (k in names(f_prod)) {
  m <- read_daily(f_prod[[k]])
  if (length(setdiff(st_names, colnames(m))))
    stop(f_prod[[k]], ": missing station columns: ",
         paste(setdiff(st_names, colnames(m)), collapse = ", "))
  mats[[k]] <- m[, st_names, drop = FALSE]
}
valid <- Reduce(`&`, lapply(mats, function(m) !is.na(m)))
mats  <- lapply(mats, function(m) { m[!valid] <- NA; m })

# ---- 3. Climatic zones -------------------------------------------------------
assign_zone <- function() {
  if (zone_method == "column") {
    if (!"Zone" %in% names(stn)) stop("No 'Zone' column in ", f_stations)
    z <- trimws(as.character(stn$Zone))
    if (!all(z %in% zone_levels)) zone_levels <<- sort(unique(z))
    return(factor(z, levels = zone_levels))
  }
  if (zone_method == "latitude")
    return(cut(stn$Latitude, c(-Inf, lat_breaks, Inf),
               labels = rev(zone_levels)) |> factor(levels = zone_levels))
  yr    <- format(all_dates, "%Y")
  tot   <- rowsum(obs_m[, st_names], yr, na.rm = TRUE)
  nval  <- rowsum(1 * !is.na(obs_m[, st_names]), yr)
  ndays <- as.vector(table(yr))
  tot[nval / ndays < min_cov_year] <- NA
  tot   <- tot / (nval / ndays)
  map   <- colMeans(tot, na.rm = TRUE)
  factor(cut(map, c(-Inf, isohyets, Inf), labels = zone_levels), levels = zone_levels)
}
zone   <- assign_zone()
n_zone <- table(zone)
print(n_zone)
zone_lab <- setNames(paste0(names(n_zone), " (n = ", as.vector(n_zone), ")"), names(n_zone))
zone_dt  <- data.table(Station = st_names,
                       Zone = factor(zone_lab[as.character(zone)], levels = zone_lab))

# ---- 4. Pentads --------------------------------------------------------------
dd  <- as.integer(format(all_dates, "%d"))
mon <- as.integer(format(all_dates, "%m"))
yy  <- as.integer(format(all_dates, "%Y"))
pid <- (yy - min(yy)) * 72L + (mon - 1L) * 6L + pmin((dd - 1L) %/% 5L + 1L, 6L)

complete <- rowsum(1 * is.na(mats$OBS), pid) == 0
comp_day <- complete[as.character(pid), , drop = FALSE]
mats     <- lapply(mats, function(m) { m[!comp_day] <- NA; m })

n_days <- length(all_dates)
dt <- rbindlist(lapply(names(mats), function(s)
  data.table(Source  = s,
             Station = rep(st_names, each = n_days),
             Day     = rep(seq_len(n_days), times = n_st),
             P       = as.vector(mats[[s]]))))
dt <- dt[!is.na(P)]
dt[, `:=`(pid = pid[Day], Month = mon[Day])]
dt[, tot := sum(P), by = .(Source, Station, pid)]

# ---- 5. Classification of the non-zero pentads -------------------------------
pm <- dt[tot > 0, {
  v  <- P[P > 0]
  nw <- length(v)
  n  <- .N
  eq <- nw >= 2L && (max(v) - min(v)) <= max(tol_abs, tol_rel * mean(v))
  .(tot = tot[1], n = n, nw = nw, Month = Month[1],
    type = if (nw == 1L) "single" else if (eq && nw == n) "full" else if (eq) "partial" else "uneven")
}, by = .(Source, Station, pid)]

pm <- zone_dt[pm, on = "Station"]
pm[, even := type %in% c("full", "partial")]
pm[, Class := cut(tot, tot_breaks, labels = tot_labels, right = TRUE)]
pm[, SourceLab := factor(src_lab[Source], levels = src_lab)]

# Console summaries
message("\nEven-split frequency (% of non-zero pentads) by series and zone:")
print(dcast(pm[, .(pct = round(100 * mean(even), 2)), by = .(SourceLab, Zone)],
            SourceLab ~ Zone, value.var = "pct"))
message("\nEven-split frequency (%) by series and month (all zones):")
print(dcast(pm[, .(pct = round(100 * mean(even), 1)), by = .(SourceLab, Month)],
            SourceLab ~ Month, value.var = "pct"))
message("\nShare of each type (% of non-zero pentads):")
print(dcast(pm[, .N, by = .(SourceLab, type)][, pct := round(100 * N / sum(N), 2), by = SourceLab],
            SourceLab ~ type, value.var = "pct"))

# ---- 6. Panel data -----------------------------------------------------------
# a) by month and zone
fa <- pm[, .(freq = 100 * mean(even), n = .N), by = .(SourceLab, Zone, Month)][n >= min_n]
fa <- fa[Month %in% plot_months]

# b) by zone and split type
fb <- pm[, .(full = 100 * mean(type == "full"), partial = 100 * mean(type == "partial")),
         by = .(SourceLab, Zone)]
fb <- melt(fb, id.vars = c("SourceLab", "Zone"), variable.name = "Type", value.name = "pct")
fb[, Type := factor(type_lab[as.character(Type)], levels = rev(type_lab))]

# c) share of the wet days of each intensity class coming from even-split pentads
dt <- pm[, .(Source, Station, pid, even)][dt, on = .(Source, Station, pid)]
dt[is.na(even), even := FALSE]
dt[, IClass := cut(P, int_breaks, labels = int_labels, right = FALSE)]
fc <- dt[P >= PMIN, .(pct = 100 * mean(even), n = .N), by = .(Source, IClass)][n >= min_n]
fc[, SourceLab := factor(src_lab[Source], levels = src_lab)]

message("\nShare (%) of wet days of each intensity class coming from even-split pentads:")
print(dcast(fc, SourceLab ~ IClass, value.var = "pct"), digits = 3)

# d) by pentad-total class
fd <- pm[, .(freq = 100 * mean(even), n = .N), by = .(SourceLab, Class)][n >= min_n]

# ---- 7. Figure ---------------------------------------------------------------
theme_fig <- theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", colour = "black", hjust = 0, size = 18),
        axis.title    = element_text(colour = "black", size = 12),
        axis.text     = element_text(colour = "black", size = 12),
        legend.title  = element_text(colour = "black", size = 12),
        legend.text   = element_text(colour = "black", size = 12),
        strip.text    = element_text(colour = "black", size = 12),
        panel.grid.minor = element_blank())

line_layers <- list(geom_line(linewidth = 0.9), geom_point(size = 2.2),
                    scale_colour_manual(values = src_cols, name = NULL))

p_a <- ggplot(fa, aes(x = Month, y = freq, colour = SourceLab, group = SourceLab)) +
  line_layers +
  facet_wrap(~ Zone, nrow = 1) +
  scale_x_continuous(breaks = plot_months, labels = substr(month.abb[plot_months], 1, 1)) +
  labs(title = "a) Even-split pentads by month",
       x = "Month", y = "Even-split pentads (% of non-zero pentads)") +
  theme_fig

p_b <- ggplot(fb, aes(x = SourceLab, y = pct, fill = Type)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  facet_wrap(~ Zone, nrow = 1) +
  scale_fill_manual(values = type_cols, name = NULL, breaks = unname(type_lab)) +
  scale_x_discrete(labels = function(x) gsub("CHIRPS ", "", x)) +
  labs(title = "b) Even-split pentads by zone and type",
       x = NULL, y = "% of non-zero pentads") +
  theme_fig

p_c <- ggplot(fc, aes(x = IClass, y = pct, colour = SourceLab, group = SourceLab)) +
  line_layers +
  labs(title = "c) Wet days generated by even splits",
       x = "Daily rainfall class (mm)", y = "% of wet days from even-split pentads") +
  theme_fig

p_d <- ggplot(fd, aes(x = Class, y = freq, colour = SourceLab, group = SourceLab)) +
  line_layers +
  labs(title = "d) Even splits by pentad total",
       x = "Pentad total (mm)", y = "Even-split pentads (%)") +
  theme_fig

fig <- (p_a / p_b / (p_c | p_d)) +
  plot_layout(guides = "collect", heights = c(1, 1, 1.1)) &
  theme(legend.position = "bottom")

ggsave(fig_file, fig, width = 16, height = 17, units = "in", dpi = 500,
       bg = "white", scale = 0.8)
message("Done: ", fig_file)