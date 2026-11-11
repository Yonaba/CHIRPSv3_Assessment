root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# 9d - Spatial organisation of daily rainfall over the 91-gauge network
#      Observations vs CHIRPS v2, CHIRPS v3 (RNL) and CHIRPS v3 (SAT), 1998-2020
# -----------------------------------------------------------------------------
# Rationale
#   Rainfall from mesoscale convective systems is patchy and organised; a
#   reanalysis-type daily driver produces widespread, spatially uniform light
#   rain. The occurrence field of each series is compared on the same days and
#   at the same stations.
#
# Method (rainy-season days only, season_months; wet = P >= occ_thr)
#   1. Occurrence correlogram: for each station pair,
#        P(both wet | one wet) = mean of P(j wet | i wet) and P(i wet | j wet)
#      against the great-circle distance. Binned means (IQR bars) and a fitted
#      decay  P(d) = c + a * exp(-d / L);  L = decorrelation length (km).
#      Dashed lines: unconditional wet-day frequency (the large-distance limit).
#   2. Daily fraction of stations reporting rain (days with at least
#      min_st_frac of the stations valid): distribution over the season, and
#      mean seasonal cycle (whole year, 15-day circular moving average).
#   SAME SAMPLE: a station-day is used only if it is valid for the 4 series.
#
# Output: figures/9d_spatial_organisation_daily_rain.png
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
fig_file   <- "figures/9d_spatial_organisation_daily_rain.png"

PMIN          <- 0.1        # mm: drizzle threshold
occ_thr       <- 1          # mm: wet-day threshold of the occurrence analysis
season_months <- 6:10        # JJAS
min_wet_days  <- 30         # min. wet days of the conditioning station for a pair
min_dist_km   <- 5          # pairs closer than this (same CHIRPS pixel) are dropped
bin_km        <- 25         # distance bin width (km)
min_pairs_bin <- 10         # min. number of pairs per distance bin
min_st_frac   <- 0.5        # min. fraction of valid stations for a daily wet fraction
smooth_days   <- 15         # moving-average window of the seasonal cycle

src_lab  <- c(OBS = "Observations", V2 = "CHIRPS v2",
              RNL = "CHIRPS v3 (RNL)", SAT = "CHIRPS v3 (SAT)")
src_cols <- c("Observations"    = "grey30",
              "CHIRPS v2"       = "darkred",
              "CHIRPS v3 (RNL)" = "darkblue",
              "CHIRPS v3 (SAT)" = "darkgreen")
season_lab <- paste(substr(month.abb[season_months], 1, 1), collapse = "")

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

mon   <- as.integer(format(all_dates, "%m"))
s_row <- mon %in% season_months

# ---- 3. Inter-station distances (km) -----------------------------------------
rad <- pi / 180
lat <- stn$Latitude * rad; lon <- stn$Longitude * rad
a   <- sin(outer(lat, lat, "-") / 2)^2 +
  outer(cos(lat), cos(lat)) * sin(outer(lon, lon, "-") / 2)^2
D   <- 2 * 6371 * asin(pmin(1, sqrt(a)))
dim(D) <- dim(a)                # pmin() drops the matrix dimensions
ut  <- which(upper.tri(D), arr.ind = TRUE)
pair_d <- D[ut]
message(sprintf("%d station pairs, distance %.0f-%.0f km", length(pair_d),
                min(pair_d), max(pair_d)))

# ---- 4. Occurrence correlogram ----------------------------------------------
cond_prob <- function(m) {
  V   <- 1 * !is.na(m)
  W   <- 1 * (!is.na(m) & m >= occ_thr)
  n11 <- crossprod(W)          # both wet
  nWV <- crossprod(W, V)       # [i, j]: i wet and j valid
  p_ji <- n11 / nWV            # P(j wet | i wet)
  p_ij <- n11 / t(nWV)         # P(i wet | j wet)
  p_ji[nWV < min_wet_days]    <- NA
  p_ij[t(nWV) < min_wet_days] <- NA
  (p_ji + p_ij) / 2
}

pairs <- rbindlist(lapply(names(mats), function(s) {
  p <- cond_prob(mats[[s]][s_row, , drop = FALSE])
  data.table(Source = s, dist = pair_d, P = p[ut])
}))[is.finite(P) & dist >= min_dist_km]
pairs[, SourceLab := factor(src_lab[Source], levels = src_lab)]

# Unconditional wet-day frequency (season)
p_wet <- rbindlist(lapply(names(mats), function(s) {
  m <- mats[[s]][s_row, , drop = FALSE]
  data.table(Source = s, p = sum(m >= occ_thr, na.rm = TRUE) / sum(!is.na(m)))
}))
p_wet[, SourceLab := factor(src_lab[Source], levels = src_lab)]

# Binned means
pairs[, bin := (floor(dist / bin_km) + 0.5) * bin_km]
binned <- pairs[, .(P = mean(P), q25 = quantile(P, 0.25), q75 = quantile(P, 0.75), n = .N),
                by = .(SourceLab, bin)][n >= min_pairs_bin]

# Exponential decay fit
fit_decay <- function(d) {
  st <- list(c = unname(quantile(d$P, 0.05)),
             a = max(0.05, unname(quantile(d$P, 0.95) - quantile(d$P, 0.05))),
             L = 100)
  tryCatch(nls(P ~ c + a * exp(-dist / L), data = d, start = st, algorithm = "port",
               lower = c(0, 0, 1), upper = c(1, 1, 5000)),
           error = function(e) NULL)
}
d_grid <- seq(0, max(pairs$dist), length.out = 300)
fits <- lapply(levels(pairs$SourceLab), function(s) {
  d <- pairs[SourceLab == s]
  f <- fit_decay(d)
  if (is.null(f)) return(list(L = NA_real_, curve = NULL))
  cf <- coef(f)
  list(L = unname(cf["L"]), c = unname(cf["c"]), a = unname(cf["a"]),
       curve = data.table(SourceLab = s, dist = d_grid,
                          P = cf["c"] + cf["a"] * exp(-d_grid / cf["L"])))
})
names(fits) <- levels(pairs$SourceLab)
curves <- rbindlist(lapply(fits, `[[`, "curve"))
if (nrow(curves)) curves[, SourceLab := factor(SourceLab, levels = src_lab)]

L_tab <- data.table(Series = names(fits),
                    L_km = round(sapply(fits, `[[`, "L")),
                    plateau_c = round(sapply(fits, function(x) if (is.null(x$c)) NA else x$c), 3),
                    wet_freq = round(p_wet$p[match(names(fits), as.character(p_wet$SourceLab))], 3))
message("\nDecorrelation lengths and wet-day frequencies (", season_lab, "):")
print(L_tab)

leg_lab <- setNames(ifelse(is.na(L_tab$L_km), L_tab$Series,
                           sprintf("%s (L = %d km)", L_tab$Series, L_tab$L_km)),
                    L_tab$Series)

# ---- 5. Daily fraction of stations reporting rain ----------------------------
wet_frac <- function(m) {
  V  <- !is.na(m)
  W  <- V & m >= occ_thr
  nv <- rowSums(V)
  f  <- rowSums(W, na.rm = TRUE) / nv
  f[nv < min_st_frac * n_st] <- NA
  f
}
fr <- rbindlist(lapply(names(mats), function(s)
  data.table(Source = s, Date = all_dates, Month = mon, f = wet_frac(mats[[s]]))))[!is.na(f)]
fr[, SourceLab := factor(src_lab[Source], levels = src_lab)]

message("\nDaily fraction of wet stations (", season_lab, "):")
print(fr[Month %in% season_months,
         .(median = round(median(f), 3), mean = round(mean(f), 3),
           pct_days_gt50 = round(100 * mean(f > 0.5), 1),
           pct_days_gt80 = round(100 * mean(f > 0.8), 1)), by = SourceLab])

# Seasonal cycle (29 February dropped, 15-day circular moving average)
md_ref <- format(seq(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "day"), "%m-%d")
fr[, md := format(Date, "%m-%d")]
cyc <- fr[md != "02-29", .(f = mean(f)), by = .(SourceLab, md)]
cyc[, doy := match(md, md_ref)]
setorder(cyc, SourceLab, doy)
smooth_circ <- function(x, k) {
  if (k <= 1) return(x)
  h  <- (k - 1) %/% 2
  xx <- c(tail(x, h), x, head(x, h))
  as.numeric(stats::filter(xx, rep(1 / k, k), sides = 2))[(h + 1):(h + length(x))]
}
cyc[, f_s := smooth_circ(f, smooth_days), by = SourceLab]

# ---- 6. Figure ---------------------------------------------------------------
theme_fig <- theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", colour = "black", hjust = 0, size = 18),
        axis.title    = element_text(colour = "black", size = 12),
        axis.text     = element_text(colour = "black", size = 12),
        legend.title  = element_text(colour = "black", size = 12),
        legend.text   = element_text(colour = "black", size = 12),
        panel.grid.minor = element_blank())

p_a <- ggplot() +
  geom_hline(data = p_wet, aes(yintercept = p, colour = SourceLab),
             linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  geom_linerange(data = binned, aes(x = bin, ymin = q25, ymax = q75, colour = SourceLab),
                 linewidth = 0.5, alpha = 0.6, show.legend = FALSE) +
  { if (nrow(curves)) geom_line(data = curves, aes(x = dist, y = P, colour = SourceLab),
                                linewidth = 1) } +
  geom_point(data = binned, aes(x = bin, y = P, colour = SourceLab), size = 2) +
  scale_colour_manual(values = src_cols, labels = leg_lab, name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.02))) +
  guides(colour = guide_legend(ncol = 2)) +
  labs(title = sprintf("a) Occurrence correlogram (%s, P \u2265 %g mm)", season_lab, occ_thr),
       x = "Inter-station distance (km)",
       y = "P(both stations wet | one wet)") +
  theme_fig +
  theme(legend.position = "bottom")

p_b <- ggplot(fr[Month %in% season_months], aes(x = f, colour = SourceLab)) +
  stat_ecdf(geom = "step", linewidth = 0.9) +
  scale_colour_manual(values = src_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(100 * x, "%")) +
  labs(title = sprintf("b) Wet-station fraction (%s)", season_lab),
       x = "Daily fraction of stations reporting rain", y = "Cumulative probability") +
  theme_fig

month_start <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
p_c <- ggplot(cyc, aes(x = doy, y = f_s, colour = SourceLab)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = src_cols, guide = "none") +
  scale_x_continuous(breaks = month_start, labels = substr(month.abb, 1, 1),
                     expand = c(0, 0)) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(title = "c) Seasonal cycle of the wet-station fraction",
       x = "Month", y = "Mean fraction of stations wet") +
  theme_fig

fig <- (p_a | (p_b / p_c)) + plot_layout(widths = c(1.25, 1))

ggsave(fig_file, fig, width = 18, height = 10, units = "in", dpi = 500,
       bg = "white", scale = 0.8)
message("Done: ", fig_file)