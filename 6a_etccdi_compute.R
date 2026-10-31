root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# =============================================================================
# 6a - ETCCDI core precipitation indices: CHIRPS v2, CHIRPS v3 (RNL), CHIRPS v3 (SAT)
# and station observations, 1998-2020.   (slow part: run once)
# -----------------------------------------------------------------------------
# Indices (climdex.pcic), as in 10_etccdi_cmp.R:
#   nP (days >= PMIN), CDD (rainy season only), CWD, SDII, R10mm, R20mm, RX1DAY,
#   RX5DAY, p90, p95, R95pTOT, R99pTOT
#
# Method
#   1. Daily CHIRPS values are extracted at the station locations.
#   2. SAME SAMPLE for the 4 series: a product day is set to missing wherever the
#      observation is missing, so differences come from the products only.
#      Values < PMIN are set to 0 (drizzle threshold), as in 10_etccdi_cmp.R.
#   3. Indices are computed per station and per year (p90 / p95: one value per
#      station, from the 1998-2020 base period). Stations are processed in
#      parallel (n_cores).
#
# Output: tables/4_etccdi_indices.csv   (long format)
#   Source, Station, year, ind, value      (year is empty for p90 / p95)
# Next step: 6b_etccdi_plot.R
# The script sets the working directory to its own folder, which must be the
# project root (the folder containing data/, gis/, figures/, tables/).
# =============================================================================

# ---- 0. Packages -------------------------------------------------------------
pkgs <- c("terra", "dplyr", "lubridate", "readr", "PCICt", "climdex.pcic")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))
library(parallel)

# ---- 1. Settings -------------------------------------------------------------
start_date <- as.Date("1998-01-01")
end_date   <- as.Date("2020-12-31")
all_dates  <- seq(start_date, end_date, by = "day")
years      <- year(start_date):year(end_date)

stations_file <- "data/csv/p_stations.csv"
obs_file      <- "data/csv/p_stations_daily_data.csv"
chirps_dir    <- "data/chirps"
out_file      <- "tables/4_etccdi_indices.csv"

products <- c("CHIRPS v2"       = "chirps_v2",
              "CHIRPS v3 (RNL)" = "chirps_v3_rnl",
              "CHIRPS v3 (SAT)" = "chirps_v3_sat")

PMIN         <- 0.1   # drizzle threshold (mm): values below are set to 0; also the nP threshold
rainy_months <- 5:10  # CDD is computed within the rainy season only (May-October)
n_cores      <- max(1, parallel::detectCores() - 1)  # 1 = no parallel computing

dir.create("tables", showWarnings = FALSE)

# ---- 2. Stations and observed daily precipitation ----------------------------
stations <- read_csv(stations_file, show_col_types = FALSE) |>
  mutate(Station = trimws(as.character(Station))) |>
  filter(!is.na(Latitude), !is.na(Longitude))

obs_wide <- read_csv(obs_file, show_col_types = FALSE, name_repair = "minimal")
names(obs_wide) <- trimws(names(obs_wide))
names(obs_wide)[1] <- "Date"
obs_dates <- as.Date(as.character(obs_wide$Date), format = "%m/%d/%Y")

common_st <- intersect(stations$Station, setdiff(names(obs_wide), "Date"))
if (length(common_st) == 0) stop("No station names match between the two csv files.")
if (length(common_st) < nrow(stations))
  message(nrow(stations) - length(common_st),
          " station(s) from p_stations.csv have no data column and are dropped.")
stations <- stations |> filter(Station %in% common_st)

# Complete daily grid (days absent from the file -> NA): days x stations
idx     <- match(all_dates, obs_dates)
obs_mat <- sapply(stations$Station, function(s)
  suppressWarnings(as.numeric(obs_wide[[s]]))[idx])
obs_mat[!is.na(obs_mat) & obs_mat < 0] <- NA  # negative values = missing-data codes

# ---- 3. CHIRPS daily series at the stations ----------------------------------
pts_ll <- vect(as.matrix(stations[, c("Longitude", "Latitude")]), crs = "EPSG:4326")

extract_chirps <- function(folder) {
  files <- list.files(file.path(chirps_dir, folder), pattern = "[0-9]{4}\\.tif$",
                      full.names = TRUE, ignore.case = TRUE)
  yrs   <- as.integer(sub(".*([0-9]{4})\\.tif$", "\\1", files, ignore.case = TRUE))
  keep  <- yrs %in% years
  files <- files[keep]; yrs <- yrs[keep]
  ord   <- order(yrs); files <- files[ord]; yrs <- yrs[ord]
  if (length(setdiff(years, yrs)) > 0)
    stop(folder, ": missing years ", paste(setdiff(years, yrs), collapse = ", "))
  
  mats <- lapply(seq_along(files), function(i) {
    r     <- rast(files[i])
    ndays <- as.integer(as.Date(paste0(yrs[i], "-12-31")) -
                          as.Date(paste0(yrs[i], "-01-01"))) + 1
    if (nlyr(r) != ndays)
      stop("File ", basename(files[i]), " has ", nlyr(r), " bands, expected ", ndays)
    pts <- if (is.lonlat(r)) pts_ll else project(pts_ll, crs(r))
    m   <- as.matrix(terra::extract(r, pts, ID = FALSE))  # stations x days
    t(m)                                                   # days x stations
  })
  m <- do.call(rbind, mats)
  if (nrow(m) != length(all_dates)) stop(folder, ": unexpected number of days")
  m[!is.na(m) & m < 0] <- NA  # NoData (-9999)
  colnames(m) <- stations$Station
  m
}

series <- list(Observations = obs_mat)
for (p in names(products)) {
  message("Extracting ", products[[p]], " ...")
  m <- extract_chirps(products[[p]])
  m[is.na(obs_mat)] <- NA  # same sample as the observations
  series[[p]] <- m
}
# Drizzle threshold
series <- lapply(series, function(m) { m[!is.na(m) & m < PMIN] <- 0; m })

# ---- 4. ETCCDI indices per station -------------------------------------------
pcict <- as.PCICt(format(all_dates, "%Y-%m-%d"), cal = "gregorian")
mon   <- as.integer(format(all_dates, "%m"))
n_yr  <- length(years)

etccdi_station <- function(x) {
  if (sum(is.finite(x)) < 365 * 3) return(NULL)  # not enough data
  
  # CDD only within the rainy season: other months are set to 1000 mm (wet)
  x_ddf <- x
  x_ddf[!(mon %in% rainy_months)] <- 1000
  
  tryCatch({
    args <- list(prec.dates = pcict, base.range = c(min(years), max(years)),
                 northern.hemisphere = TRUE, prec.qtiles = c(0.9, 0.95, 0.99))
    ci     <- do.call(climdexInput.raw, c(list(prec = x),     args))
    ci_ddf <- do.call(climdexInput.raw, c(list(prec = x_ddf), args))
    
    ann <- list(
      nP      = climdex.rnnmm(ci, PMIN),
      CDD     = climdex.cdd(ci_ddf, spells.can.span.years = FALSE),
      CWD     = climdex.cwd(ci, spells.can.span.years = FALSE),
      SDII    = climdex.sdii(ci),
      R10mm   = climdex.r10mm(ci),
      R20mm   = climdex.r20mm(ci),
      RX1DAY  = climdex.rx1day(ci, freq = "annual"),
      RX5DAY  = climdex.rx5day(ci, freq = "annual"),
      R95pTOT = climdex.r95ptot(ci),
      R99pTOT = climdex.r99ptot(ci)
    )
    vals <- lapply(ann, as.numeric)
    if (any(lengths(vals) != n_yr)) stop("unexpected number of annual values")
    
    list(annual = data.frame(year  = rep(years, times = length(vals)),
                             ind   = rep(names(vals), each = n_yr),
                             value = unlist(vals, use.names = FALSE)),
         static = c(p90 = ci@quantiles[["prec"]][["q90"]],
                    p95 = ci@quantiles[["prec"]][["q95"]]))
  }, error = function(e) NULL)
}

# One task per (series, station)
tasks <- expand.grid(Source = names(series), Station = stations$Station,
                     stringsAsFactors = FALSE)
xs <- lapply(seq_len(nrow(tasks)), function(i)
  series[[tasks$Source[i]]][, tasks$Station[i]])

message(nrow(tasks), " series to process on ", n_cores, " core(s) ...")
if (n_cores > 1) {
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("etccdi_station", "pcict", "mon", "years", "n_yr",
                      "PMIN", "rainy_months"))
  invisible(clusterEvalQ(cl, { library(PCICt); library(climdex.pcic) }))
  outs <- tryCatch(
    parLapplyLB(cl, xs, function(x) suppressWarnings(etccdi_station(x))),
    finally = stopCluster(cl))
} else {
  outs <- lapply(xs, function(x) suppressWarnings(etccdi_station(x)))
}

n_failed <- sum(vapply(outs, is.null, logical(1)))
if (n_failed > 0)
  message(n_failed, " series skipped (not enough data or computation error).")

# ---- 5. Save -------------------------------------------------------------------
res <- bind_rows(lapply(seq_along(outs), function(i) {
  out <- outs[[i]]
  if (is.null(out)) return(NULL)
  bind_rows(
    data.frame(Source = tasks$Source[i], Station = tasks$Station[i],
               year = out$annual$year, ind = out$annual$ind, value = out$annual$value),
    data.frame(Source = tasks$Source[i], Station = tasks$Station[i],
               year = NA_integer_, ind = names(out$static), value = unname(out$static))
  )
}))

write_csv(res, out_file)
message("Done: ", out_file, " (", nrow(res), " rows)")