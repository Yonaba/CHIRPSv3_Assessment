Sys.setenv(TZ = "UTC")

library(terra)
library(SPEI)
library(ggplot2)
library(ggpubr)

# ---------------------------------------------------------------- settings
y.start <- 1998
y.end   <- 2020
spi.scales <- c(1, 3, 6, 12)

products <- c("Observations", "CHIRPS v2", "CHIRPS v3 (RNL)", "CHIRPS v3 (SAT)")
chirps.products <- products[-1]
chirps.dirs <- c("CHIRPS v2"        = "data/chirps/chirps_v2",
                 "CHIRPS v3 (RNL)"  = "data/chirps/chirps_v3_rnl",
                 "CHIRPS v3 (SAT)"  = "data/chirps/chirps_v3_sat")

# colour convention
cols <- c("Observations"   = "black",
          "CHIRPS v2"       = "darkred",
          "CHIRPS v3 (RNL)" = "darkblue",
          "CHIRPS v3 (SAT)" = "darkgreen")
ltys <- c("Observations"   = "dotted",
          "CHIRPS v2"       = "solid",
          "CHIRPS v3 (RNL)" = "solid",
          "CHIRPS v3 (SAT)" = "solid")

dir.create("figures", showWarnings = FALSE)

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

# ---------------------------------------------------------------- stations
stations <- read.csv("data/csv/p_stations.csv", header = TRUE, stringsAsFactors = FALSE)
obs.raw  <- read.csv("data/csv/p_stations_daily_data.csv", header = TRUE,
                     check.names = FALSE, stringsAsFactors = FALSE)

obs.dates <- as.Date(obs.raw[[1]], format = "%m/%d/%Y")
obs.mat   <- as.matrix(obs.raw[, -1, drop = FALSE])
colnames(obs.mat) <- colnames(obs.raw)[-1]

# keep stations present in both tables, in the same order
st.names <- intersect(stations$Station, colnames(obs.mat))
if (length(st.names) == 0) stop("No matching station names between the two csv files.")
stations <- stations[match(st.names, stations$Station), ]
obs.mat  <- obs.mat[, st.names, drop = FALSE]

# ---------------------------------------------------------------- daily date axis (1998-2020)
all.dates <- seq(as.Date(paste0(y.start, "-01-01")), as.Date(paste0(y.end, "-12-31")), by = "day")

# observations re-indexed on the full daily axis (missing days -> NA)
obs.daily <- obs.mat[match(all.dates, obs.dates), , drop = FALSE]

# ---------------------------------------------------------------- CHIRPS extraction at stations
pts.ll <- vect(data.frame(lon = stations$Longitude, lat = stations$Latitude),
               geom = c("lon", "lat"), crs = "EPSG:4326")

extract.chirps <- function(dir) {
  files <- list.files(dir, pattern = "[0-9]{4}\\.tif$", full.names = TRUE)
  yrs <- as.integer(sub(".*([0-9]{4})\\.tif$", "\\1", files))
  keep <- yrs >= y.start & yrs <= y.end
  files <- files[keep]; yrs <- yrs[keep]
  files <- files[order(yrs)]; yrs <- sort(yrs)
  if (!all(y.start:y.end %in% yrs)) stop("Missing years in ", dir)
  
  out <- vector("list", length(files))
  for (i in seq_along(files)) {
    r <- rast(files[i])
    ndays <- as.integer(as.Date(paste0(yrs[i], "-12-31")) - as.Date(paste0(yrs[i], "-01-01"))) + 1
    if (nlyr(r) != ndays) stop(basename(files[i]), ": ", nlyr(r), " bands, expected ", ndays)
    pts <- if (same.crs(crs(r), "EPSG:4326")) pts.ll else project(pts.ll, crs(r))
    v <- as.matrix(terra::extract(r, pts, ID = FALSE))     # stations x days
    v[v < 0] <- NA                                  # CHIRPS no-data (-9999)
    out[[i]] <- t(v)                                # days x stations
  }
  m <- do.call(rbind, out)
  colnames(m) <- stations$Station
  m
}

chirps.daily <- lapply(chirps.dirs, function(d) {
  message("Extracting: ", d)
  extract.chirps(d)
})

# ---------------------------------------------------------------- monthly precipitation (station-average)
ym <- format(all.dates, "%Y-%m")
month.dates <- as.Date(paste0(unique(ym), "-01"))

monthly.avg <- function(daily) {
  # monthly total per station (NA if any day is missing), then mean across stations
  ms <- rowsum(daily, group = ym, reorder = FALSE)
  m <- rowMeans(ms, na.rm = TRUE)
  m[is.nan(m)] <- NA
  as.numeric(m)
}

pr <- data.frame(Date = month.dates,
                 "Observations"    = monthly.avg(obs.daily),
                 "CHIRPS v2"       = monthly.avg(chirps.daily[["CHIRPS v2"]]),
                 "CHIRPS v3 (RNL)" = monthly.avg(chirps.daily[["CHIRPS v3 (RNL)"]]),
                 "CHIRPS v3 (SAT)" = monthly.avg(chirps.daily[["CHIRPS v3 (SAT)"]]),
                 check.names = FALSE)

# ---------------------------------------------------------------- SPI
get.spi <- function(x, s) {
  xt <- ts(x, start = c(y.start, 1), frequency = 12)
  v <- as.numeric(spi(xt, s, na.rm = TRUE, verbose = FALSE)$fitted)
  v[!is.finite(v)] <- NA
  v
}

spi.list <- lapply(spi.scales, function(s) {
  sapply(products, function(p) get.spi(pr[[p]], s))
})
names(spi.list) <- spi.scales

# ---------------------------------------------------------------- correlation helper
cor.fun <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(c(rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value)
}
fmt.p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", fmt(p, 3)))
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

# ---------------------------------------------------------------- Figure 1: SPI time series
base.theme <- theme_bw(base_size = 12) +
  theme(plot.title   = element_text(face = "bold", hjust = 0, colour = "black", size = 18),
        axis.text    = element_text(size = 12, colour = "black"),
        axis.title   = element_text(size = 12, colour = "black"),
        legend.text  = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"))

ts.panels <- list()
for (i in seq_along(spi.scales)) {
  s <- spi.scales[i]
  m <- spi.list[[as.character(s)]]
  
  long <- do.call(rbind, lapply(products, function(p)
    data.frame(Date = pr$Date, Product = p, SPI = m[, p])))
  long$Product <- factor(long$Product, levels = products)
  
  # correlation annotations (CHIRPS vs Observations)
  ann <- data.frame(Product = chirps.products, y = c(3.3, 2.9, 2.5))
  ann$label <- sapply(chirps.products, function(p) {
    r <- cor.fun(m[, p], m[, "Observations"])
    # plotmath expression (parsed by annotate) so that the rho symbol is displayed
    paste0("'", p, ":'~rho~'= ", fmt(r["rho"]), ", ", fmt.p(r["p"]), "'")
  })
  
  g <- ggplot(long, aes(x = Date, y = SPI, colour = Product, linetype = Product)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(data = subset(long, Product != "Observations"), linewidth = 0.6, na.rm = TRUE) +
    geom_line(data = subset(long, Product == "Observations"), linewidth = 0.8, na.rm = TRUE) +
    scale_colour_manual(values = cols, breaks = products) +
    scale_linetype_manual(values = ltys, breaks = products) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y",
                 limits = as.Date(c(paste0(y.start, "-01-01"), paste0(y.end, "-12-31"))),
                 expand = expansion(mult = 0.01)) +
    coord_cartesian(ylim = c(-3.5, 3.5)) +
    labs(title = paste0(letters[i], ") SPI ", s, "-month", ifelse(s == 1, "", "s")),
         x = NULL, y = if (i %in% c(1, 3)) "SPI [-]" else NULL) +   # y title: left-most panels only
    base.theme +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # x-axis text: bottom-most panels only
  if (!(i %in% c(3, 4))) {
    g <- g + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  for (k in seq_len(nrow(ann))) {
    g <- g + annotate("text", x = as.Date(paste0(y.start, "-03-01")), y = ann$y[k],
                      label = ann$label[k], colour = cols[[ann$Product[k]]],
                      hjust = 0, size = 4, parse = TRUE)
  }
  ts.panels[[i]] <- g
}

fig1 <- ggarrange(plotlist = ts.panels, nrow = 2, ncol = 2, align = "hv",
                  common.legend = TRUE, legend = "bottom")

ggsave("figures/7_spi_timeseries.png", plot = fig1,
       width = 20, height = 12, units = "in", dpi = 500, scale = 0.8, bg = "white")

# ---------------------------------------------------------------- Figure 2: correlation heatmaps
res <- do.call(rbind, lapply(spi.scales, function(s) {
  m <- spi.list[[as.character(s)]]
  do.call(rbind, lapply(chirps.products, function(p) {
    x <- m[, p]; o <- m[, "Observations"]
    full <- cor.fun(x, o)
    dry  <- cor.fun(ifelse(o < 0, x, NA), ifelse(o < 0, o, NA))   # months where observed SPI < 0
    wet  <- cor.fun(ifelse(o > 0, x, NA), ifelse(o > 0, o, NA))   # months where observed SPI > 0
    data.frame(Scale = s, Product = p,
               Condition = c("Full series", "Dry condition (SPI < 0)", "Wet condition (SPI > 0)"),
               rho = c(full["rho"], dry["rho"], wet["rho"]),
               p   = c(full["p"],   dry["p"],   wet["p"]))
  }))
}))
res$Product <- factor(res$Product, levels = chirps.products)
res$ScaleF  <- factor(paste0(res$Scale, "-month", ifelse(res$Scale == 1, "", "s")),
                      levels = rev(paste0(spi.scales, "-month", ifelse(spi.scales == 1, "", "s"))))
res$label   <- ifelse(is.na(res$rho), "", paste0(fmt(res$rho), stars(res$p)))

conds  <- c("Full series", "Dry condition (SPI < 0)", "Wet condition (SPI > 0)")
titles <- c("Complete SPI series", "Dry condition (SPI < 0)", "Wet condition (SPI > 0)")
rho.lim <- c(0.6, 1)

hm.panels <- lapply(seq_along(conds), function(i) {
  d <- subset(res, Condition == conds[i])
  ggplot(d, aes(x = Product, y = ScaleF, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = label), colour = "black", size = 4.5) +
    # three-colour gradient (values below 0.6 are squished to the lowest colour)
    scale_fill_gradientn(colours = c("#ffffb2", "#fd8d3c", "#e31a1c"),
                         limits = rho.lim, breaks = seq(0.6, 1, 0.1),
                         oob = scales::squish, na.value = "grey85",
                         name = expression("Spearman's"~rho)) +
    guides(fill = guide_colourbar(barwidth = unit(18, "lines"), barheight = unit(1, "lines"))) +
    scale_x_discrete(expand = c(0, 0),
                     labels = function(x) sub(" \\(", "\n(", x)) +   # wrap "CHIRPS v3 (RNL)" on two lines
    scale_y_discrete(expand = c(0, 0)) +
    labs(title = paste0(letters[i], ") ", titles[i]), x = NULL,
         y = if (i == 1) "Timescale" else NULL) +   # y title: leftmost panel only
    base.theme +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5))
})

fig2 <- ggarrange(plotlist = hm.panels, nrow = 1, ncol = 3, align = "hv",
                  common.legend = TRUE, legend = "bottom")

ggsave("figures/7_spi_correlation_heatmaps.png", plot = fig2,
       width = 22, height = 7, units = "in", dpi = 500, scale = 0.8, bg = "white")

print("finished.")