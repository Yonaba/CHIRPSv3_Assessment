root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# 9a - Within-pentad partitioning of daily rainfall
#      Observations vs CHIRPS v2, CHIRPS v3 (RNL) and CHIRPS v3 (SAT), 1998-2020
# -----------------------------------------------------------------------------
# Rationale
#   CHIRPS v3 RNL and SAT share the same pentadal totals: they only differ by the
#   daily "driver" used to split each pentad into days (reanalysis for RNL,
#   satellite for SAT). The within-pentad partitioning is therefore the daily
#   signature of the driver, and can be diagnosed from the station series alone.
#
# Method
#   1. Daily series at the 91 stations (CSV files written by
#      1_extract_stations_daily.R). SAME SAMPLE for the 4 series: a day is kept
#      only if it is valid for all of them. Values < PMIN are set to 0.
#   2. CHIRPS pentads: 6 per month (days 1-5, 6-10, 11-15, 16-20, 21-25 and
#      26-end of month, i.e. 3 to 6 days). Only complete pentads are used.
#   3. For each station-pentad with a non-zero total, and for each series:
#        CR   = share of the pentad total falling on the wettest day
#        NW   = number of wet days (P >= PMIN)
#        Gini = normalised Gini index of the daily shares (printed only)
#               (0 = evenly spread over all days, 1 = all on one day)
#        H    = normalised Shannon entropy (printed only)
#      stratified by pentad-total class (each series' own total) and zone.
#   4. Isolated-storm composite: gauge days with P >= storm_min mm carrying at
#      least storm_cr of the gauge rainfall of the surrounding 7-day window;
#      mean rainfall of each series at lags -3..+3 days.
#   5. One example pentad (single gauge storm, largest RNL - SAT wet-day gap).
#
# Climatic zones (zone_method):
#   "isohyet"  : observed mean annual rainfall (< 600 Sahelian, 600-900
#                Sudano-Sahelian, > 900 Sudanian)
#   "latitude" : > 14 N Sahelian, 11.5-14 N Sudano-Sahelian, < 11.5 N Sudanian
#   "column"   : a "Zone" column of p_stations.csv
#
# Output: figures/9a_within_pentad_partitioning.png
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
fig_file   <- "figures/9a_within_pentad_partitioning.png"

PMIN         <- 0.1             # mm: drizzle threshold, also the wet-day threshold
zone_method  <- "isohyet"       # "isohyet", "latitude" or "column"
isohyets     <- c(600, 900)     # mm/year
lat_breaks   <- c(11.5, 14)     # degrees N
min_cov_year <- 0.90            # min. fraction of valid days for a year to enter the annual mean

tot_breaks <- c(0, 5, 15, 30, 60, Inf)              # pentad-total classes (mm)
tot_labels <- c("0-5", "5-15", "15-30", "30-60", ">60")

storm_min  <- 20                # mm: gauge storm day
storm_cr   <- 0.70              # min. share of the 7-day gauge total falling on the storm day
lags       <- -3:3
ex_min_tot <- 30                # example pentad: min. gauge pentad total (mm)
ex_min_cr  <- 0.85              # example pentad: min. gauge concentration ratio

src_lab  <- c(OBS = "Observations", V2 = "CHIRPS v2",
              RNL = "CHIRPS v3 (RNL)", SAT = "CHIRPS v3 (SAT)")
src_cols <- c("Observations"    = "grey45",
              "CHIRPS v2"       = "darkred",
              "CHIRPS v3 (RNL)" = "darkblue",
              "CHIRPS v3 (SAT)" = "darkgreen")
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
  m[!is.na(m) & m < 0]    <- NA   # missing-data codes
  m[!is.na(m) & m < PMIN] <- 0    # drizzle threshold
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

# Same sample for the 4 series
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
  tot   <- tot / (nval / ndays)        # scale slightly incomplete years to 365 days
  map   <- colMeans(tot, na.rm = TRUE)
  factor(cut(map, c(-Inf, isohyets, Inf), labels = zone_levels), levels = zone_levels)
}
zone   <- assign_zone()
n_zone <- table(zone)
print(n_zone)
zone_lab <- setNames(paste0(names(n_zone), " (n = ", as.vector(n_zone), ")"), names(n_zone))
zone_dt  <- data.table(Station = st_names,
                       Zone = factor(zone_lab[as.character(zone)], levels = zone_lab))

# ---- 4. CHIRPS pentads -------------------------------------------------------
dd  <- as.integer(format(all_dates, "%d"))
mon <- as.integer(format(all_dates, "%m"))
yy  <- as.integer(format(all_dates, "%Y"))
pid <- (yy - min(yy)) * 72L + (mon - 1L) * 6L + pmin((dd - 1L) %/% 5L + 1L, 6L)

# Complete pentads only (same for the 4 series, since the sample is shared)
complete <- rowsum(1 * is.na(mats$OBS), pid) == 0          # pentads x stations
comp_day <- complete[as.character(pid), , drop = FALSE]     # days x stations
mats_d   <- mats                                            # daily same-sample series
mats_p   <- lapply(mats, function(m) { m[!comp_day] <- NA; m })

# Check: RNL and SAT pentad totals
T_rnl <- rowsum(mats_p$RNL, pid); T_sat <- rowsum(mats_p$SAT, pid)
dT    <- abs(T_rnl - T_sat)
message(sprintf(paste0("RNL vs SAT pentad totals (after the %g mm drizzle threshold): ",
                       "%.2f%% of pentads agree within max(0.1 mm, 1%%) (max diff = %.2f mm)"),
                PMIN, 100 * mean(dT <= pmax(0.1, 0.01 * T_rnl), na.rm = TRUE),
                max(dT, na.rm = TRUE)))

# ---- 5. Partitioning metrics per station-pentad ------------------------------
n_days <- length(all_dates)
dt <- rbindlist(lapply(names(mats_p), function(s)
  data.table(Source  = s,
             Station = rep(st_names, each = n_days),
             Day     = rep(seq_len(n_days), times = n_st),
             P       = as.vector(mats_p[[s]]))))
dt <- dt[!is.na(P)]
dt[, `:=`(pid = pid[Day], Month = mon[Day])]
dt[, tot := sum(P), by = .(Source, Station, pid)]

pm <- dt[tot > 0, {
  n <- .N
  s <- sort(P / tot)
  .(tot   = tot[1],
    n     = n,
    Month = Month[1],
    CR    = s[n],
    NW    = sum(P >= PMIN),
    Gini  = sum((2 * seq_len(n) - n - 1) * s) / (n - 1),
    H     = -sum(s[s > 0] * log(s[s > 0])) / log(n))
}, by = .(Source, Station, pid)]

pm <- zone_dt[pm, on = "Station"]
pm[, Class := cut(tot, tot_breaks, labels = tot_labels, right = TRUE)]
pm[, SourceLab := factor(src_lab[Source], levels = src_lab)]

# Console summary
message("\nMedian partitioning metrics by series and zone:")
print(pm[, .(n_pentads = .N, mean_NW = round(mean(NW), 2),
             median_CR = round(median(CR), 3), median_Gini = round(median(Gini), 3),
             median_H = round(median(H), 3),
             pct_all_days_wet = round(100 * mean(NW == n), 1)),
         by = .(SourceLab, Zone)][order(Zone, SourceLab)])

# ---- 6. Isolated-storm composite --------------------------------------------
shift_m <- function(m, k) {
  out <- matrix(NA_real_, nrow(m), ncol(m))
  if (k >= 0) out[seq_len(n_days - k), ] <- m[(1 + k):n_days, ]
  else        out[(1 - k):n_days, ]      <- m[seq_len(n_days + k), ]
  out
}
wsum  <- Reduce(`+`, lapply(lags, function(k) shift_m(mats_d$OBS, k)))
storm <- !is.na(wsum) & !is.na(mats_d$OBS) & mats_d$OBS >= storm_min &
  mats_d$OBS / wsum >= storm_cr
storm[is.na(storm)] <- FALSE
n_storm <- sum(storm)
message("\nIsolated gauge storms in the composite: ", n_storm)

comp <- rbindlist(lapply(names(mats_d), function(s) {
  rbindlist(lapply(lags, function(k) {
    x <- shift_m(mats_d[[s]], k)[storm]
    data.table(Source = s, lag = k, mean = mean(x, na.rm = TRUE),
               se = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  }))
}))
comp[, SourceLab := factor(src_lab[Source], levels = src_lab)]
print(dcast(comp, lag ~ SourceLab, value.var = "mean"), digits = 3)

# ---- 7. Example pentad -------------------------------------------------------
w <- dcast(pm, Station + pid ~ Source, value.var = c("tot", "CR", "NW"))
cand <- w[!is.na(tot_OBS) & tot_OBS >= ex_min_tot & CR_OBS >= ex_min_cr &
            !is.na(NW_RNL) & !is.na(NW_SAT)]
cand[, score := NW_RNL - NW_SAT]
ex <- if (nrow(cand)) cand[order(-score, -tot_OBS)][1] else NULL

# ---- 8. Figure ---------------------------------------------------------------
theme_fig <- theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", colour = "black", hjust = 0, size = 18),
        plot.subtitle = element_text(colour = "black", size = 12, hjust = 0),
        axis.title    = element_text(colour = "black", size = 12),
        axis.text     = element_text(colour = "black", size = 12),
        legend.title  = element_text(colour = "black", size = 12),
        legend.text   = element_text(colour = "black", size = 12),
        strip.text    = element_text(colour = "black", size = 12),
        panel.grid.minor = element_blank())

box_panel <- function(var, ttl, ylab) {
  ggplot(pm, aes(x = Class, y = .data[[var]], fill = SourceLab)) +
    geom_boxplot(outlier.shape = NA, width = 0.75, linewidth = 0.35, alpha = 0.85,
                 position = position_dodge(0.85)) +
    stat_summary(aes(group = SourceLab), fun = mean, geom = "point", shape = 23,
                 size = 1.8, fill = "white", colour = "black",
                 position = position_dodge(0.85)) +
    facet_wrap(~ Zone, nrow = 1) +
    scale_fill_manual(values = src_cols, name = NULL) +
    labs(title = ttl, x = "Pentad total (mm)", y = ylab) +
    theme_fig
}

p_a <- box_panel("CR", "a) Share of the pentad total on the wettest day",
                 "Concentration ratio (-)") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
p_b <- box_panel("NW", "b) Number of wet days per pentad",
                 paste0("Wet days (P \u2265 ", PMIN, " mm)")) +
  scale_y_continuous(breaks = 0:6)

p_c <- ggplot(comp, aes(x = lag, y = mean, colour = SourceLab)) +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                width = 0.15, linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = src_cols, name = NULL) +
  scale_x_continuous(breaks = lags) +
  guides(colour = "none") +
  labs(title = "c) Isolated-storm composite",
       subtitle = sprintf("n = %s gauge storms (P \u2265 %g mm, \u2265 %g%% of the 7-day total)",
                          format(n_storm, big.mark = ","), storm_min, 100 * storm_cr),
       x = "Days from the gauge storm day", y = "Mean rainfall (mm/day)") +
  theme_fig

if (!is.null(ex)) {
  idx <- which(pid == ex$pid)
  ex_d <- rbindlist(lapply(names(mats_d), function(s)
    data.table(Source = s, Date = all_dates[idx], P = mats_d[[s]][idx, ex$Station])))
  ex_d[, SourceLab := factor(src_lab[Source], levels = src_lab)]
  ex_d[, DayLab := factor(format(Date, "%d %b"), levels = format(all_dates[idx], "%d %b"))]
  p_d <- ggplot(ex_d, aes(x = DayLab, y = P, fill = SourceLab)) +
    geom_col(position = position_dodge(0.85), width = 0.8, colour = "black",
             linewidth = 0.2) +
    scale_fill_manual(values = src_cols, name = NULL) +
    guides(fill = "none") +
    labs(title = "d) Example of a single gauge storm",
         subtitle = sprintf("%s, %s | wet days: Obs = %d, v2 = %d, RNL = %d, SAT = %d",
                            ex$Station, format(all_dates[idx[1]], "%Y"),
                            ex$NW_OBS, ex$NW_V2, ex$NW_RNL, ex$NW_SAT),
         x = NULL, y = "Daily rainfall (mm)") +
    theme_fig
} else {
  p_d <- ggplot() + annotate("text", x = 0, y = 0, label = "No pentad matches the criteria") +
    labs(title = "d) Example of a single gauge storm") + theme_void() +
    theme(plot.title = element_text(face = "bold", colour = "black", hjust = 0, size = 18))
}

fig <- (p_a / p_b / (p_c | p_d)) +
  plot_layout(guides = "collect", heights = c(1, 1, 1.15)) &
  theme(legend.position = "bottom")

ggsave(fig_file, fig, width = 16, height = 16.5, units = "in", dpi = 500,
       bg = "white", scale = 0.8)
message("Done: ", fig_file)