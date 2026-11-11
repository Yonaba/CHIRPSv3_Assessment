root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# 9c - Disaggregation counterfactuals: pentad-total error vs daily-timing error
#      CHIRPS v2, CHIRPS v3 (RNL) and CHIRPS v3 (SAT), 1998-2020
# -----------------------------------------------------------------------------
# Each CHIRPS daily value = pentad total x daily ratio (driver). Swapping the
# ratios isolates the two sources of daily error.
#
# Series (same daily sample, complete CHIRPS pentads only)
#   CHIRPS v2, CHIRPS v3 (RNL), CHIRPS v3 (SAT)   actual products
#   v2 totals + gauge timing   v2 pentad totals x gauge daily ratios
#   v3 totals + gauge timing   v3 pentad totals x gauge daily ratios
#                              ("oracle": perfect daily timing; the remaining
#                               daily error comes from the pentad totals only)
#   v3 totals + v2 timing      v3 pentad totals x CHIRPS v2 daily ratios
#                              (cross-disaggregation)
#   When the donor pentad total is 0 (no ratio available), the pentad is split
#   evenly (1/n), which is the least-squares split without timing information.
#
# Decomposition of the daily mean squared error (MSE) of each product
#   MSE(product) = MSE(oracle)                    -> "Pentad totals"
#                + [MSE(product) - MSE(oracle)]   -> "Daily timing (driver)"
#   and, for the cross-disaggregation, MSE(v3 totals + v2 timing) - MSE(v3 oracle)
#   = timing error of the v2 driver applied to the v3 totals.
#   NB: the oracle keeps the gauge timing, not the least-squares split, so the
#   timing term is a "perfect-timing" difference and may be slightly negative
#   at a few stations.
#
# Station-level scores (all valid days of eval_months): KGE (Gupta et al., 2009),
# Pearson r, RMSE, critical success index (CSI) for P >= occ_thr.
#
# Output: figures/9c_disaggregation_counterfactuals.png
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
fig_file   <- "figures/9c_disaggregation_counterfactuals.png"

PMIN         <- 0.1         # mm: drizzle threshold
occ_thr      <- 1           # mm: rain/no-rain threshold of the CSI
eval_months  <- 6:10        # months used for the station scores (e.g. 5:10 = rainy season)
plot_months  <- 6:10        # months shown in panel d
min_days     <- 365 * 3     # min. number of valid days per station
zone_method  <- "isohyet"   # "isohyet", "latitude" or "column" (see 9a)
isohyets     <- c(600, 900)
lat_breaks   <- c(11.5, 14)
min_cov_year <- 0.90
same_tot_tol <- 0.99        # RNL and SAT totals considered identical if >= 99% of
# pentads agree within max(0.1 mm, 1% of the total)

zone_levels <- c("Sahelian", "Sudano-Sahelian", "Sudanian")

dir.create("figures", showWarnings = FALSE)

# ---- 2. Data -----------------------------------------------------------------
# Product values are read raw (no drizzle threshold) so that the pentad totals,
# and hence the counterfactuals, are exactly those of the products.
read_daily <- function(f, zero_drizzle = TRUE) {
  if (!file.exists(f)) stop(f, " not found (run 1_extract_stations_daily.R first).")
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  names(d) <- trimws(names(d)); names(d)[1] <- "Date"
  dts <- as.Date(as.character(d$Date), format = "%m/%d/%Y")
  x <- d[match(all_dates, dts), -1, drop = FALSE]
  m <- matrix(suppressWarnings(as.numeric(unlist(x, use.names = FALSE))),
              nrow = nrow(x), dimnames = list(NULL, names(x)))
  m[!is.na(m) & m < 0]    <- NA
  if (zero_drizzle) m[!is.na(m) & m < PMIN] <- 0
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
  m <- read_daily(f_prod[[k]], zero_drizzle = FALSE)
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
zone_dt  <- data.table(Station = st_names, Zone = zone_lab[as.character(zone)])

# ---- 4. Pentad totals, daily ratios and counterfactual series ----------------
dd  <- as.integer(format(all_dates, "%d"))
mon <- as.integer(format(all_dates, "%m"))
yy  <- as.integer(format(all_dates, "%Y"))
pid <- (yy - min(yy)) * 72L + (mon - 1L) * 6L + pmin((dd - 1L) %/% 5L + 1L, 6L)
pid_c  <- as.character(pid)
n_days <- length(all_dates)

# Pentad totals expanded back to days (NA if the pentad is incomplete)
pent_tot <- lapply(mats, function(m) rowsum(m, pid))                 # pentads x stations
day_tot  <- lapply(pent_tot, function(p) p[pid_c, , drop = FALSE])   # days x stations
n_in_pent <- matrix(as.vector(table(pid))[match(pid, sort(unique(pid)))],
                    nrow = n_days, ncol = n_st)

dT <- abs(pent_tot$RNL - pent_tot$SAT)
share_same <- mean(dT <= pmax(0.1, 0.01 * pent_tot$RNL), na.rm = TRUE)
message(sprintf("RNL vs SAT pentad totals: %.2f%% of pentads agree within max(0.1 mm, 1%%) (max diff = %.2f mm)",
                100 * share_same, max(dT, na.rm = TRUE)))
same_tot <- share_same >= same_tot_tol
if (!same_tot)
  warning("RNL and SAT pentad totals differ: counterfactuals are shown for each variant.")

daily_ratio <- function(m, Tm) {
  r <- m / Tm
  z <- !is.na(Tm) & Tm == 0
  r[z] <- 1 / n_in_pent[z]            # no donor timing: even split
  r
}
r_obs <- daily_ratio(mats$OBS, day_tot$OBS)
r_v2  <- daily_ratio(mats$V2,  day_tot$V2)

incomplete <- is.na(day_tot$OBS)      # same NA pattern for every series
obs <- mats$OBS; obs[incomplete] <- NA

S <- list(
  V2     = mats$V2,
  V2_OR  = day_tot$V2  * r_obs,
  RNL    = mats$RNL,
  SAT    = mats$SAT,
  RNL_OR = day_tot$RNL * r_obs,
  RNL_X  = day_tot$RNL * r_v2,
  SAT_OR = day_tot$SAT * r_obs,
  SAT_X  = day_tot$SAT * r_v2
)
S <- lapply(S, function(m) { m[incomplete] <- NA; colnames(m) <- st_names; m })

ser_lab <- if (same_tot) {
  c(V2 = "CHIRPS v2", V2_OR = "v2 totals + gauge timing",
    RNL = "CHIRPS v3 (RNL)", SAT = "CHIRPS v3 (SAT)",
    RNL_OR = "v3 totals + gauge timing", RNL_X = "v3 totals + v2 timing")
} else {
  c(V2 = "CHIRPS v2", V2_OR = "v2 totals + gauge timing",
    RNL = "CHIRPS v3 (RNL)", RNL_OR = "RNL totals + gauge timing",
    RNL_X = "RNL totals + v2 timing",
    SAT = "CHIRPS v3 (SAT)", SAT_OR = "SAT totals + gauge timing",
    SAT_X = "SAT totals + v2 timing")
}
ser_cols <- c(V2 = "darkred", V2_OR = "#e7a4a4", RNL = "darkblue", SAT = "darkgreen",
              RNL_OR = "orange3", RNL_X = "purple4", SAT_OR = "gold3", SAT_X = "orchid4")
ser_cols <- setNames(ser_cols[names(ser_lab)], ser_lab)

# ---- 5. Station-level scores -------------------------------------------------
keep_rows <- mon %in% eval_months

station_scores <- function(sim) {
  rbindlist(lapply(seq_len(n_st), function(j) {
    s <- sim[keep_rows, j]; o <- obs[keep_rows, j]
    ok <- !is.na(s) & !is.na(o)
    if (sum(ok) < min_days) return(NULL)
    s <- s[ok]; o <- o[ok]
    r  <- suppressWarnings(cor(s, o))
    be <- mean(s) / mean(o)
    ga <- (sd(s) / mean(s)) / (sd(o) / mean(o))   # CV ratio (Kling et al., 2012)
    h  <- sum(s >= occ_thr & o >= occ_thr)
    f  <- sum(s >= occ_thr & o <  occ_thr)
    mi <- sum(s <  occ_thr & o >= occ_thr)
    data.table(Station = st_names[j], n = sum(ok), r = r,
               KGE  = 1 - sqrt((r - 1)^2 + (ga - 1)^2 + (be - 1)^2),
               MSE  = mean((s - o)^2), RMSE = sqrt(mean((s - o)^2)),
               CSI  = h / (h + f + mi))
  }))
}

sc <- rbindlist(lapply(names(S), function(k) station_scores(S[[k]])[, Series := k]))
sc <- zone_dt[sc, on = "Station"]

message("\nMedian station scores:")
print(sc[Series %in% names(ser_lab),
         .(n_stations = .N, KGE = round(median(KGE), 3), r = round(median(r), 3),
           RMSE = round(median(RMSE), 2), CSI = round(median(CSI), 3)),
         by = Series][, Series := ser_lab[Series]][])

# ---- 6. MSE decomposition ----------------------------------------------------
decomp_def <- data.table(
  Product = c("CHIRPS v2", "CHIRPS v3 (RNL)", "CHIRPS v3 (SAT)",
              if (same_tot) "v3 totals + v2 timing" else c("RNL totals + v2 timing", "SAT totals + v2 timing")),
  act     = c("V2", "RNL", "SAT", if (same_tot) "RNL_X" else c("RNL_X", "SAT_X")),
  orc     = c("V2_OR", "RNL_OR", "SAT_OR", if (same_tot) "RNL_OR" else c("RNL_OR", "SAT_OR")))
prod_levels <- decomp_def$Product

mse_z <- rbind(sc[, .(MSE = mean(MSE)), by = .(Series, Zone)],
               sc[, .(Zone = "All stations", MSE = mean(MSE)), by = Series],
               use.names = TRUE)
dec <- rbindlist(lapply(seq_len(nrow(decomp_def)), function(i) {
  a <- mse_z[Series == decomp_def$act[i]]
  o <- mse_z[Series == decomp_def$orc[i]]
  m <- merge(a, o, by = "Zone", suffixes = c("_act", "_orc"))
  m[, .(Zone, Product = decomp_def$Product[i],
        `Pentad totals` = MSE_orc, `Daily timing (driver)` = MSE_act - MSE_orc,
        total = MSE_act)]
}))
dec[, timing_share := 100 * `Daily timing (driver)` / total]

message("\nDaily MSE decomposition (mm2/day2):")
print(dec[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

dec_long <- melt(dec, id.vars = c("Zone", "Product", "total", "timing_share"),
                 variable.name = "Component", value.name = "MSE")
dec_long[, Component := factor(Component, levels = c("Daily timing (driver)", "Pentad totals"))]
dec_long[, Product := factor(Product, levels = prod_levels)]
dec_long[, Zone := factor(Zone, levels = c(unname(zone_lab), "All stations"))]
dec_lab <- unique(dec_long[, .(Zone, Product, total, timing_share)])

# ---- 7. Timing share of the MSE by month (all stations pooled) ---------------
mon_mat <- matrix(mon, nrow = n_days, ncol = n_st)
mse_month <- function(sim) {
  e  <- (sim - obs)^2
  ok <- !is.na(e)
  data.table(Month = mon_mat[ok], e = e[ok])[, .(mse = mean(e)), by = Month]
}
mm <- rbindlist(lapply(unique(c(decomp_def$act, decomp_def$orc)),
                       function(k) mse_month(S[[k]])[, Series := k]))
dm <- rbindlist(lapply(seq_len(nrow(decomp_def)), function(i) {
  a <- mm[Series == decomp_def$act[i]]; o <- mm[Series == decomp_def$orc[i]]
  m <- merge(a, o, by = "Month", suffixes = c("_act", "_orc"))
  m[, .(Month, Product = decomp_def$Product[i],
        share = 100 * (mse_act - mse_orc) / mse_act)]
}))
dm <- dm[Month %in% plot_months]
dm[, Product := factor(Product, levels = prod_levels)]

prod_cols <- c("CHIRPS v2" = "darkred", "CHIRPS v3 (RNL)" = "darkblue",
               "CHIRPS v3 (SAT)" = "darkgreen", "v3 totals + v2 timing" = "purple4",
               "RNL totals + v2 timing" = "purple4", "SAT totals + v2 timing" = "orchid4")
prod_cols <- prod_cols[prod_levels]

# ---- 8. Figure ---------------------------------------------------------------
theme_fig <- theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", colour = "black", hjust = 0, size = 18),
        axis.title    = element_text(colour = "black", size = 12),
        axis.text     = element_text(colour = "black", size = 12),
        legend.title  = element_text(colour = "black", size = 12),
        legend.text   = element_text(colour = "black", size = 12),
        strip.text    = element_text(colour = "black", size = 12),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

sc_plot <- sc[Series %in% names(ser_lab)]
sc_plot[, SerLab := factor(ser_lab[Series], levels = ser_lab)]

score_panel <- function(var, ttl, ylab) {
  ggplot(sc_plot, aes(x = SerLab, y = .data[[var]], fill = SerLab)) +
    geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.35, alpha = 0.85) +
    geom_jitter(width = 0.12, height = 0, size = 0.8, alpha = 0.35, colour = "black") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2.2,
                 fill = "white", colour = "black") +
    scale_fill_manual(values = ser_cols, guide = "none") +
    labs(title = ttl, x = NULL, y = ylab) +
    theme_fig +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
}

p_a <- score_panel("KGE", "a) Daily KGE", "KGE (-)")
p_b <- score_panel("CSI", paste0("b) Daily CSI (P \u2265 ", occ_thr, " mm)"), "CSI (-)")

p_c <- ggplot(dec_long, aes(x = Product, y = MSE, fill = Component)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_text(data = dec_lab, inherit.aes = FALSE,
            aes(x = Product, y = total, label = sprintf("%.0f%%", timing_share)),
            vjust = -0.4, size = 3.6, colour = "black") +
  facet_wrap(~ Zone, nrow = 1) +
  scale_fill_manual(values = c("Daily timing (driver)" = "#d95f02",
                               "Pentad totals" = "grey65"), name = NULL) +
  scale_x_discrete(labels = function(x) gsub("CHIRPS ", "", x)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "c) Decomposition of the daily MSE (label: timing share)",
       x = NULL, y = expression(MSE~(mm^2~day^{-2}))) +
  theme_fig +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

p_d <- ggplot(dm, aes(x = Month, y = share, colour = Product, group = Product)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = prod_cols, name = NULL) +
  scale_x_continuous(breaks = plot_months, labels = month.abb[plot_months]) +
  labs(title = "d) Timing share of the daily MSE by month",
       x = "Month", y = "Daily timing share of MSE (%)") +
  theme_fig

fig <- ((p_a | p_b) / p_c / p_d) + plot_layout(heights = c(1.1, 1, 0.9))

ggsave(fig_file, fig, width = 16, height = 19, units = "in", dpi = 500,
       bg = "white", scale = 0.8)
message("Done: ", fig_file)