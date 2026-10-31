root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# Wavelet analysis of observed and CHIRPS rainfall (biwavelet package)
# -----------------------------------------------------------------------------
# Same plotting logic as 8_biwavelet.R (continuous wavelet transform, normalised
# power, significance contours, cone of influence), adapted to this case.
#
# Panels
#   - Wavelet power spectra: Observations, CHIRPS v2, CHIRPS v3 (RNL), CHIRPS v3 (SAT)
#   - Wavelet coherence (with phase arrows): observations vs each product
#   - Global wavelet spectrum (mean power per period) of the 4 series
#
# Combined figure (2 x 4 + a colour-bar column), columns aligned by dataset:
#     (a) Global spectrum | (b) coherence v2 | (c) coherence RNL | (d) coherence SAT | bar
#     (e) Obs. power      | (f) power v2     | (g) power RNL     | (h) power SAT      | bar
#   (rows are swapped with coherence_on_top <- FALSE)
#
# Alignment: the colour bars are NOT drawn by biwavelet (plot.cb = FALSE, which
# shrinks the plotting area of the panel that carries the bar). They are drawn
# by draw_cb() in a dedicated, narrow layout column, so all the panels of a
# column have exactly the same plotting area in both rows.
#
# Series: regional mean monthly rainfall (mm/month), 1998-2020.
#   1. CHIRPS daily values are extracted at the station locations.
#   2. A COMMON mask is applied: a (station, day) is used only if the observation
#      AND the three products are all valid, so the 4 series share the same
#      sampling and differences between panels come from the products only.
#   3. Monthly totals per station (month kept if >= min_cov_month of its days valid),
#      then the mean over stations (month kept if >= min_stations stations).
#   4. Remaining empty months (identical for the 4 series) are linearly
#      interpolated, since the wavelet transform needs a regular series.
#   5. Analysis series: square root of the monthly totals (reduces the skew from
#      wet-season peaks), then removal of the mean seasonal cycle (calendar-month
#      climatology of each series). The annual cycle no longer hides the other bands.
#
# Output (figures/): 5_biwavelet_combined.png
#   (the three separate figures are only written if save_individual = TRUE)
# The script sets the working directory to its own folder, which must be the
# project root (the folder containing data/, gis/, figures/, tables/).
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(biwavelet)
  library(viridis)
})

# ---- 0. Settings -------------------------------------------------------------
f_stations <- "data/csv/p_stations.csv"
f_obs      <- "data/csv/p_stations_daily_data.csv"
dir_chirps <- "data/chirps"
dir_fig    <- "figures"
fig_combined  <- "5_biwavelet_combined.png"
fig_power     <- "5_biwavelet_power.png"
fig_coherence <- "5_biwavelet_coherence.png"
fig_global    <- "5_biwavelet_global_spectrum.png"
save_individual  <- FALSE  # TRUE = also save the three separate figures
coherence_on_top <- TRUE   # combined figure: coherence row first (FALSE = power row first)

year_start <- 1998
year_end   <- 2020

min_cov_month <- 0.90  # min. fraction of valid days for a station-month to be kept
min_stations  <- 10    # min. number of stations for a regional-mean month to be kept

# Preprocessing of the regional monthly series
use_sqrt        <- TRUE  # square-root transform of monthly totals
remove_seasonal <- TRUE  # subtract the mean seasonal cycle (calendar-month means)

# Wavelet settings (periods in months)
siglevel     <- 0.95
minp         <- 2
maxp         <- 32
lperiods     <- c(2, 4, 8, 16, 32)  # axis labels
zlim_power   <- c(-3, 3)   # colour scale of the power panels, log2(power / variance);
# NULL = derived from the data. Values outside are clamped.
zlim_coh     <- c(0.4, 1)  # colour scale of the coherence panels (values below are clamped)
alpha_coi    <- 0.15       # opacity of the cone-of-influence shading (default 0.5)
lwd_sig      <- 1.2        # line width of the significance contours
nrands_wtc   <- 300        # Monte Carlo runs for the coherence significance (lower = faster)
arrow_cutoff <- 0.9        # phase arrows are drawn where coherence >= this value
min_valid_frac <- 1/3      # global spectrum: min. fraction of times outside the cone of influence

# Text sizes (all increased; the margins below were enlarged accordingly)
cex_axis   <- 2.3   # tick labels (also used for the colour bars)
cex_lab    <- 2.5   # axis titles
cex_main   <- 2.6   # panel titles
cex_legend <- 1.8   # legend / annotation of the global spectrum (top-right corner)
mgp_fig    <- c(5, 1.2, 0)  # distance of axis title / tick labels from the axis (lines)

# Combined figure: width / height of one plotting area (> 1 = wider than tall).
# Lower it to make the figure taller, raise it to make it flatter.
panel_aspect <- 1.3

products <- data.frame(
  id     = c("chirps_v2", "chirps_v3_rnl", "chirps_v3_sat"),
  label  = c("CHIRPS v2", "CHIRPS v3 (RNL)", "CHIRPS v3 (SAT)"),
  colour = c("#8B0000", "#00008B", "#006400"),  # dark red, dark blue, dark green
  stringsAsFactors = FALSE
)

dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
set.seed(42)  # the significance test is Monte Carlo

days_in_year <- function(y) ifelse((y %% 4 == 0 & y %% 100 != 0) | y %% 400 == 0, 366, 365)

# ---- 1. Observations (wide -> long) -----------------------------------------
stations <- read.csv(f_stations, check.names = FALSE, strip.white = TRUE,
                     stringsAsFactors = FALSE)
stations$Station <- trimws(stations$Station)

obs_wide <- read.csv(f_obs, check.names = FALSE, strip.white = TRUE,
                     stringsAsFactors = FALSE)
names(obs_wide) <- trimws(names(obs_wide))
names(obs_wide)[1] <- "Date"
obs_wide$Date <- as.Date(obs_wide$Date, format = "%m/%d/%Y")
yr <- as.integer(format(obs_wide$Date, "%Y"))
obs_wide <- obs_wide[!is.na(yr) & yr >= year_start & yr <= year_end, ]

obs_stations <- setdiff(names(obs_wide), "Date")
obs <- tidyr::pivot_longer(obs_wide, -Date, names_to = "Station", values_to = "OBS") |>
  mutate(OBS = suppressWarnings(as.numeric(OBS)),
         OBS = ifelse(OBS < 0, NA, OBS))  # negative values = missing-data codes

st <- stations |>
  filter(Station %in% obs_stations, !is.na(Latitude), !is.na(Longitude))
message(nrow(st), " stations with coordinates and daily data.")

# ---- 2. CHIRPS files and extraction at the stations --------------------------
index_files <- function(p) {
  f <- list.files(file.path(dir_chirps, p), pattern = "[0-9]{4}\\.tif$",
                  full.names = TRUE, ignore.case = TRUE)
  d <- data.frame(file = f,
                  year = as.integer(sub(".*([0-9]{4})\\.tif$", "\\1", f, ignore.case = TRUE)))
  d <- d[d$year >= year_start & d$year <= year_end, ]
  d[order(d$year), ]
}
chirps_files <- setNames(lapply(products$id, index_files), products$id)
for (p in products$id) {
  if (nrow(chirps_files[[p]]) == 0) stop("No tif file found for ", p)
  miss <- setdiff(year_start:year_end, chirps_files[[p]]$year)
  if (length(miss)) warning(p, ": missing years: ", paste(miss, collapse = ", "))
}

xy         <- as.matrix(st[, c("Longitude", "Latitude")])
pts_chirps <- vect(xy, crs = "EPSG:4326")  # station points (same order as st)

extract_chirps <- function(fdf) {
  dplyr::bind_rows(lapply(seq_len(nrow(fdf)), function(i) {
    r <- rast(fdf$file[i])
    y <- fdf$year[i]
    n <- nlyr(r)
    if (n != days_in_year(y))
      warning(basename(fdf$file[i]), ": ", n, " bands (expected ", days_in_year(y), ")")
    # Dates from the band names (e.g. "19980101"); fallback: consecutive days from 1 Jan
    dates <- as.Date(names(r), format = "%Y%m%d")
    if (anyNA(dates) || length(dates) != n)
      dates <- seq(as.Date(sprintf("%d-01-01", y)), by = "day", length.out = n)
    m <- as.matrix(terra::extract(r, pts_chirps, ID = FALSE))  # stations x days
    if (ncol(m) != n || nrow(m) != nrow(st))
      stop(basename(fdf$file[i]), ": extracted ", nrow(m), " x ", ncol(m),
           " values, expected ", nrow(st), " x ", n)
    m[m < 0] <- NA  # no-data codes (e.g. -9999)
    data.frame(Station = rep(st$Station, times = n),
               Date    = rep(dates, each = nrow(st)),
               value   = as.vector(m))
  }))
}

daily <- obs |> filter(Station %in% st$Station)
for (p in products$id) {
  message("Extracting ", p, " ...")
  cl <- extract_chirps(chirps_files[[p]])
  names(cl)[names(cl) == "value"] <- p
  daily <- left_join(daily, cl, by = c("Station", "Date"))
}

# ---- 3. Common mask, monthly totals, regional mean ---------------------------
vars  <- c("OBS", products$id)
daily <- as.data.frame(daily)
ok    <- complete.cases(daily[, vars])
daily[!ok, vars] <- NA
daily$ok    <- ok
daily$month <- as.Date(format(daily$Date, "%Y-%m-01"))

monthly <- daily |>
  group_by(Station, month) |>
  summarise(n_ok = sum(ok),
            across(all_of(vars), \(x) sum(x, na.rm = TRUE)),
            .groups = "drop") |>
  mutate(n_days = as.integer(as.Date(format(month + 32, "%Y-%m-01")) - month)) |>
  filter(n_ok / n_days >= min_cov_month)

areal <- monthly |>
  group_by(month) |>
  summarise(n_st = n(), across(all_of(vars), mean), .groups = "drop") |>
  filter(n_st >= min_stations)

months_all <- seq(as.Date(sprintf("%d-01-01", year_start)),
                  as.Date(sprintf("%d-12-01", year_end)), by = "month")
ts_df <- left_join(data.frame(month = months_all), areal, by = "month")
n     <- nrow(ts_df)
t_idx <- seq_len(n)  # time axis in months (dt = 1)

miss <- is.na(ts_df$OBS)
message(sum(miss), " of ", n, " months have no regional mean and are linearly interpolated.")
if (mean(miss) > 0.10)
  warning("More than 10% of the months are interpolated: consider relaxing ",
          "min_cov_month / min_stations.")
for (v in vars) {
  y <- ts_df[[v]]
  if (any(miss)) y[miss] <- approx(t_idx[!miss], y[!miss], xout = t_idx[miss], rule = 2)$y
  ts_df[[v]] <- y
}

# ---- 3b. Analysis series: sqrt transform + removal of the seasonal cycle -----
mon <- as.integer(format(ts_df$month, "%m"))
ana <- ts_df
for (v in vars) {
  y <- ts_df[[v]]
  if (use_sqrt) y <- sqrt(y)
  if (remove_seasonal) y <- y - ave(y, mon, FUN = mean)  # minus calendar-month mean
  ana[[v]] <- y
}

# ---- 4. Wavelet transforms and coherence ----------------------------------------
panels <- data.frame(
  id     = vars,
  label  = c("Observations", products$label),
  colour = c("black", products$colour),
  stringsAsFactors = FALSE
)

CWTs <- lapply(seq_len(nrow(panels)), function(k) {
  message("Wavelet transform: ", panels$label[k])
  wt(cbind(t = t_idx, pr = ana[[panels$id[k]]]), sig.level = siglevel, max.scale = maxp)
})
names(CWTs) <- panels$id

WTC <- lapply(seq_len(nrow(products)), function(j) {
  message("Wavelet coherence: Observations vs ", products$label[j])
  wtc(cbind(t = t_idx, pr = ana$OBS),
      cbind(t = t_idx, pr = ana[[products$id[j]]]),
      sig.level = siglevel, max.scale = maxp, nrands = nrands_wtc)
})

# TRUE where a wavelet coefficient lies outside the cone of influence
valid_mask <- function(w) {
  if (length(w$coi) == ncol(w$power)) {
    outer(w$period, w$coi, "<=")
  } else {
    warning("Unexpected length of the cone of influence: it is not excluded.")
    matrix(TRUE, nrow(w$power), ncol(w$power))
  }
}

# Variance used to normalise the power (as in type = "power.norm")
power_var <- function(w, y) if (!is.null(w$sigma2)) w$sigma2 else var(y)

# Colour scale of the power panels, common to the 4 panels
if (is.null(zlim_power)) {
  zv <- unlist(lapply(seq_along(CWTs), function(k) {
    z <- log2(CWTs[[k]]$power / power_var(CWTs[[k]], ana[[panels$id[k]]]))
    z[valid_mask(CWTs[[k]])]
  }))
  zv         <- zv[is.finite(zv)]
  zlim_power <- c(floor(unname(quantile(zv, 0.01))), ceiling(unname(quantile(zv, 0.995))))
}
message("Power colour scale (log2 normalised power): ", zlim_power[1], " to ", zlim_power[2])

# Global wavelet spectrum: mean power per period over time, outside the cone of
# influence; periods with less than min_valid_frac of valid times are not drawn.
global_spec <- function(w) {
  m  <- valid_mask(w)
  pw <- w$power
  pw[!m] <- NA
  g  <- rowMeans(pw, na.rm = TRUE)
  g[rowMeans(m) < min_valid_frac] <- NA
  g
}
G   <- sapply(CWTs, global_spec)   # periods x series
per <- CWTs[[1]]$period

# ---- 5. Panel drawing functions (draw into the current device) -----------------
pal      <- turbo(64, begin = 0.2)
yr_all   <- year_start:year_end
x_at     <- seq(1, n, by = 24)                 # January of every second year
x_labels <- yr_all[seq(1, length(yr_all), by = 2)]
lwd_g    <- c(2.8, 2, 2, 2)                    # global spectrum: obs, v2, RNL, SAT
lty_g    <- c(1, 1, 1, 2)                      # SAT dashed so it stays visible over RNL
ylab_per <- "Period (months)"

add_axes <- function() {
  axis(side = 1, at = x_at, labels = x_labels, las = 2, cex.axis = cex_axis)
  axis(side = 2, at = log2(lperiods), labels = lperiods, cex.axis = cex_axis)
}

# NB: tol = 1 draws the significance contour exactly at `siglevel`; the
# significance threshold itself is set by sig.level in wt() / wtc().

# show_cb = TRUE lets biwavelet draw its own colour bar (only used by the
# individual figures). In the combined figure it is always FALSE and the bars
# come from draw_cb(), so that all the plotting areas are identical.
draw_power <- function(k, mar, show_ylab, show_cb, title_txt) {
  w <- CWTs[[k]]
  # Clamp to the colour scale so that no cell is left blank (white)
  s2 <- power_var(w, ana[[panels$id[k]]])
  w$power <- pmin(pmax(w$power, s2 * 2^zlim_power[1]), s2 * 2^zlim_power[2])
  
  par(mar = mar)
  plot(w, tol = 1, ylim = c(minp, maxp), type = "power.norm",
       plot.cb = show_cb, plot.phase = FALSE,
       fill.cols = pal, zlim = zlim_power,
       lty.sig = 1, lwd.sig = lwd_sig,
       lty.coi = 1, lwd.coi = 2, col.coi = "black", alpha.coi = alpha_coi,
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = if (show_ylab) ylab_per else "", cex.lab = cex_lab)
  title(title_txt, line = 1.0, adj = 0, cex.main = cex_main, col.main = panels$colour[k])
  add_axes()
}

draw_coh <- function(j, mar, show_ylab, show_cb, title_txt) {
  W <- WTC[[j]]
  W$rsq <- pmin(pmax(W$rsq, zlim_coh[1]), zlim_coh[2])  # clamp to the colour scale
  
  par(mar = mar)
  plot(W, tol = 1, ylim = c(minp, maxp), plot.phase = TRUE, arrow.cutoff = arrow_cutoff,
       plot.cb = show_cb, fill.cols = pal, zlim = zlim_coh,
       lty.sig = 1, lwd.sig = lwd_sig,
       lty.coi = 1, lwd.coi = 2, col.coi = "black", alpha.coi = alpha_coi,
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = if (show_ylab) ylab_per else "", cex.lab = cex_lab)
  title(title_txt, line = 1.0, adj = 0, cex.main = cex_main, col.main = products$colour[j])
  add_axes()
}

# Stand-alone vertical colour bar, drawn in its own layout cell. Same top and
# bottom margins as the panels, so it spans exactly the height of the plots.
draw_cb <- function(zl, at, labels, mar) {
  par(mar = mar)
  zz   <- seq(zl[1], zl[2], length.out = length(pal) + 1)
  zmid <- (zz[-1] + zz[-length(zz)]) / 2
  image(x = c(0, 1), y = zz, z = matrix(zmid, nrow = 1), zlim = zl, col = pal,
        xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  axis(side = 4, at = at, labels = labels, las = 1, cex.axis = cex_axis)
  box()
}
draw_cb_power <- function(mar)
  draw_cb(zlim_power, at = seq(zlim_power[1], zlim_power[2]),
          labels = format(2^seq(zlim_power[1], zlim_power[2]), trim = TRUE, digits = 2),
          mar = mar)
draw_cb_coh <- function(mar)
  draw_cb(zlim_coh, at = seq(zlim_coh[1], zlim_coh[2], by = 0.1),
          labels = format(seq(zlim_coh[1], zlim_coh[2], by = 0.1), nsmall = 1),
          mar = mar)

# Global spectrum with the period on the y-axis, so that it shares the axis of
# the wavelet panels next to it
draw_global <- function(mar, show_ylab, title_txt) {
  par(mar = mar)
  plot(NA, xlim = c(0, max(G, na.rm = TRUE) * 2.2), ylim = log2(c(maxp, minp)),
       yaxs = "i", yaxt = "n", xlab = "Global wavelet power",
       ylab = if (show_ylab) ylab_per else "", cex.lab = cex_lab, cex.axis = cex_axis)
  axis(side = 2, at = log2(lperiods), labels = lperiods, cex.axis = cex_axis)
  for (k in seq_len(nrow(panels)))
    lines(G[, k], log2(per), col = panels$colour[k], lwd = lwd_g[k], lty = lty_g[k])
  legend("topright", legend = panels$label, col = panels$colour, lwd = lwd_g, lty = lty_g,
         bty = "n", cex = cex_legend, text.font = 2, seg.len = 1.6)
  title(title_txt, line = 1.0, adj = 0, cex.main = cex_main)
}

open_png <- function(file, w, h) {
  png(file.path(dir_fig, file), units = "in", width = w, height = h, res = 600, pointsize = 8)
  par(font.axis = 2, font.lab = 2, cex.axis = cex_axis, mgp = mgp_fig)  # cex.axis also reaches the colour bars
}

# ---- 6. Combined figure: 2 rows x (4 panels + 1 colour-bar column) -------------
# Margins (in lines) per column
mL <- c(8.5, 5, 5, 5)
mR <- rep(1.5, 4)
mb <- 7.5   # bottom margin (rotated year labels)
mt <- 5     # top margin (panel titles)
mar_c <- function(c) c(mb, mL[c], mt, mR[c])
tag   <- function(i) paste0("(", letters[i], ") ")

# Colour-bar column: left gap, bar width, right margin (tick labels)
cb_gap <- 1; cb_bar <- 2.2; cb_lab <- 6.5
mar_cb <- c(mb, cb_gap, mt, cb_lab)

# Figure size: the width is fixed; the height follows from panel_aspect so that
# each plotting area has the requested width / height ratio.
line_in <- 8 * 1.2 / 72                       # inches per margin line at pointsize 8
fig_w   <- 20
inner_w <- (fig_w / line_in - sum(mL + mR) - (cb_gap + cb_bar + cb_lab)) /
  length(mL) * line_in                        # one plotting area, inches
fig_h   <- 2 * (inner_w / panel_aspect + (mb + mt) * line_in)
message(sprintf("Combined figure: %.1f x %.1f in", fig_w, fig_h))

open_png(fig_combined, fig_w, fig_h)
# Column widths chosen so that the four plotting areas have equal widths;
# 5th column = colour bars
layout(matrix(1:10, nrow = 2, byrow = TRUE),
       widths = c((inner_w / line_in) + mL + mR, cb_gap + cb_bar + cb_lab))
par(font.axis = 2, font.lab = 2, cex.axis = cex_axis, mgp = mgp_fig)

row_coherence <- function(t0) {
  draw_global(mar_c(1), TRUE, paste0(tag(t0), "Global spectrum"))
  for (j in seq_len(nrow(products)))
    draw_coh(j, mar_c(j + 1), FALSE, FALSE,
             paste0(tag(t0 + j), "Obs. vs ", products$label[j]))
  draw_cb_coh(mar_cb)
}
row_power <- function(t0) {
  for (k in seq_len(nrow(panels)))
    draw_power(k, mar_c(k), k == 1, FALSE,
               paste0(tag(t0 + k - 1), panels$label[k]))
  draw_cb_power(mar_cb)
}

if (coherence_on_top) {
  row_coherence(1); row_power(5)
} else {
  row_power(1); row_coherence(5)
}
dev.off()

# ---- 7. Individual figures ------------------------------------------------------
if (save_individual) {
  # Power spectra, 2 x 2
  open_png(fig_power, 10, 9)
  par(mfrow = c(2, 2))
  for (k in seq_len(nrow(panels)))
    draw_power(k, c(8, ifelse(k %in% c(1, 3), 9, 5), 5, ifelse(k %in% c(2, 4), 12, 2)),
               k %in% c(1, 3), k %in% c(2, 4), paste0(tag(k), panels$label[k]))
  dev.off()
  
  # Coherence, 1 x 3
  open_png(fig_coherence, 15, 6)
  par(mfrow = c(1, 3))
  for (j in seq_len(nrow(products)))
    draw_coh(j, c(8, ifelse(j == 1, 9, 5), 5, ifelse(j == 3, 12, 2)),
             j == 1, j == 3, paste0(tag(j), "Observations vs ", products$label[j]))
  dev.off()
  
  # Global spectrum
  open_png(fig_global, 5.5, 6.5)
  draw_global(c(8, 9, 5, 2), TRUE, "Global wavelet spectrum")
  dev.off()
}

message("Done: figures written to ", dir_fig)