root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# ------------------------------------------------------------------
# Extract CHIRPS daily precipitation at station locations (1998-2020)
# and export one CSV per product, with the same structure as
# data/csv/p_stations_daily_data.csv
# ------------------------------------------------------------------

library(terra)

# ---- Settings ----------------------------------------------------
start_date <- as.Date("1998-01-01")
end_date   <- as.Date("2020-12-31")
years      <- 1998:2020

stations_file <- "data/csv/p_stations.csv"
obs_file      <- "data/csv/p_stations_daily_data.csv"
chirps_dir    <- "data/chirps"
out_dir       <- "data/csv"

products <- c(
  chirpsv2    = "chirps_v2",
  chirpsv3rnl = "chirps_v3_rnl",
  chirpsv3sat = "chirps_v3_sat"
)

# ---- Load station coordinates -----------------------------------
stations <- read.csv(stations_file, stringsAsFactors = FALSE, check.names = FALSE)
stations$Station <- trimws(as.character(stations$Station))

# ---- Load observed data (only to reproduce the exact structure) --
obs <- read.csv(obs_file, stringsAsFactors = FALSE, check.names = FALSE)
obs_names <- names(obs)[-1]            # station names, in the original order
obs_names <- trimws(obs_names)

# Keep stations that are in both files, ordered as in the observed file
missing_in_coords <- setdiff(obs_names, stations$Station)
if (length(missing_in_coords) > 0) {
  warning("Stations without coordinates (will be filled with NA): ",
          paste(missing_in_coords, collapse = ", "))
}
stations <- stations[match(obs_names, stations$Station), ]
stations$Station <- obs_names

# Station points (WGS84, decimal degrees)
has_xy <- !is.na(stations$Longitude) & !is.na(stations$Latitude)
pts <- vect(stations[has_xy, c("Longitude", "Latitude")],
            geom = c("Longitude", "Latitude"), crs = "EPSG:4326")

# Complete daily date sequence
all_dates <- seq(start_date, end_date, by = "day")

# ---- Process each product ---------------------------------------
for (prod in names(products)) {
  
  message("Processing ", prod, " ...")
  folder <- file.path(chirps_dir, products[[prod]])
  
  files <- list.files(folder, pattern = "[0-9]{4}\\.tif$",
                      full.names = TRUE, ignore.case = TRUE)
  file_years <- as.integer(sub(".*([0-9]{4})\\.tif$", "\\1", files,
                               ignore.case = TRUE))
  files <- files[order(file_years)]
  file_years <- sort(file_years)
  keep <- file_years %in% years
  files <- files[keep]
  file_years <- file_years[keep]
  
  year_list <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    yr <- file_years[i]
    r  <- rast(files[i])
    
    # Station coordinates reprojected if the raster is not in lon/lat
    p <- if (!is.na(crs(r)) && !is.lonlat(r)) project(pts, crs(r)) else pts
    
    vals <- terra::extract(r, p, ID = FALSE)    # n_stations x n_bands
    vals <- as.matrix(vals)
    
    # Dates for this year's bands
    yr_dates <- seq(as.Date(sprintf("%d-01-01", yr)),
                    as.Date(sprintf("%d-12-31", yr)), by = "day")
    if (nlyr(r) != length(yr_dates)) {
      stop(sprintf("%s: %d bands but %d days in %d",
                   basename(files[i]), nlyr(r), length(yr_dates), yr))
    }
    
    m <- t(vals)                                # days x stations
    rownames(m) <- format(yr_dates, "%Y-%m-%d")
    year_list[[i]] <- m
    rm(r, vals, m); gc(verbose = FALSE)
    message("  ", yr, " done")
  }
  
  sim <- do.call(rbind, year_list)
  
  # CHIRPS no-data (-9999) -> NA
  sim[sim < 0] <- NA
  
  # Full matrix on the complete date sequence, all stations of the obs file
  out_mat <- matrix(NA_real_, nrow = length(all_dates), ncol = length(obs_names),
                    dimnames = list(format(all_dates, "%Y-%m-%d"), obs_names))
  idx_cols <- which(has_xy)
  rows_ok  <- intersect(rownames(sim), rownames(out_mat))
  out_mat[rows_ok, idx_cols] <- sim[rows_ok, , drop = FALSE]
  
  # Same structure as the observed file: Date ("m/d/Y") + station columns
  out <- data.frame(Date = format(all_dates, "%m/%d/%Y"),
                    out_mat, check.names = FALSE, stringsAsFactors = FALSE)
  # Match the observed date format exactly (no zero-padding, if that is how
  # the observed file is written)
  if (!any(grepl("^0|/0", obs[[1]][1:min(50, nrow(obs))]))) {
    out$Date <- paste0(as.integer(format(all_dates, "%m")), "/",
                       as.integer(format(all_dates, "%d")), "/",
                       format(all_dates, "%Y"))
  }
  names(out)[1] <- names(obs)[1]
  names(out)[-1] <- names(obs)[-1]
  
  out_file <- file.path(out_dir, paste0("p_", prod, "_daily_data.csv"))
  write.csv(out, out_file, row.names = FALSE, na = "NA")
  message("  -> written: ", out_file)
}

message("All done.")