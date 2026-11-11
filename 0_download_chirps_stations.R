root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_folder)
Sys.setenv(TZ = "UTC")

# -----------------------------------------------------------------------------
# Download CHIRPS monthly "stations used" diagnostic files (1998-2020)
# for CHIRPS v3 and CHIRPS v2, then build one aggregated CSV per version,
# filtered to Burkina Faso.
#
#   v3: https://data.chc.ucsb.edu/products/CHIRPS/v3.0/diagnostics/
#       list_of_stations_used/monthly/global.stationsUsed.v3.YYYY.MM.csv
#   v2: https://data.chc.ucsb.edu/products/CHIRPS/v2.0/diagnostics/
#       list_of_stations_used/monthly/global.stationsUsed.YYYY.MM.csv
#       + extra.stationsUsed.YYYY.MM.csv (stations added after the 2015
#         release; only published for 2016.01 - 2019.12)
# -----------------------------------------------------------------------------

library(data.table)

# ---- Settings ---------------------------------------------------------------
years          <- 1998:2020
months         <- 1:12
stations_dir   <- file.path("data", "chirps_stations")
target_country <- "Burkina_Faso"
max_tries      <- 3       # retries per file
overwrite      <- FALSE   # re-download files already on disk?

# Each version lists its file patterns. Optional per-type date windows
# (YYYYMM, inclusive) restrict which months are requested from the server;
# types without a window are requested for every year x month.
versions <- list(
  v3 = list(
    base_url = "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/diagnostics/list_of_stations_used/monthly",
    patterns = c(global = "global.stationsUsed.v3.%d.%02d.csv"),
    windows  = list(),
    optional = character(0)
  ),
  v2 = list(
    base_url = "https://data.chc.ucsb.edu/products/CHIRPS/v2.0/diagnostics/list_of_stations_used/monthly",
    patterns = c(global = "global.stationsUsed.%d.%02d.csv",
                 extra  = "extra.stationsUsed.%d.%02d.csv"),
    windows  = list(extra = c(from = 201601L, to = 201912L)),
    optional = "extra"    # tolerate an occasional gap inside the window
  )
)

options(timeout = max(500, getOption("timeout")))

# ---- Helper: download one file with retries ---------------------------------
download_one <- function(url, dest, max_tries = 3) {
  tmp <- paste0(dest, ".part")
  for (i in seq_len(max_tries)) {
    ok <- tryCatch({
      status <- download.file(url, tmp, mode = "wb", quiet = TRUE)
      status == 0 && file.exists(tmp) && file.size(tmp) > 0
    }, error = function(e) FALSE, warning = function(w) FALSE)
    
    if (isTRUE(ok)) {
      file.rename(tmp, dest)  # only keep complete downloads
      return(TRUE)
    }
    if (file.exists(tmp)) unlink(tmp)
    Sys.sleep(2 * i)          # back off before retrying
  }
  FALSE
}

# ---- Helper: keep only the target country -----------------------------------
# Matches "Burkina_Faso" and "Burkina Faso" alike, in case versions differ.
filter_country <- function(dt, country) {
  col <- grep("^country_name$", names(dt), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(col)) {
    warning("No 'country_name' column found; file skipped.")
    return(NULL)
  }
  norm <- function(x) tolower(gsub("[ _]+", "_", trimws(x)))
  dt[norm(get(col)) == norm(country)]
}

# ---- Helper: build the list of files to request -----------------------------
build_grid <- function(cfg, out_dir) {
  grid <- CJ(year = years, month = months, type = names(cfg$patterns))
  grid[, ym := year * 100L + month]
  
  # Drop months outside each type's publication window (no server request)
  for (tp in names(cfg$windows)) {
    w <- cfg$windows[[tp]]
    grid <- grid[type != tp | (ym >= w[["from"]] & ym <= w[["to"]])]
  }
  
  grid[, fname := sprintf(cfg$patterns[type], year, month), by = type]
  grid[, url   := file.path(cfg$base_url, fname)]
  grid[, dest  := file.path(out_dir, fname)]
  grid[, status := NA_character_]
  grid[, ym := NULL]
  setorder(grid, year, month, type)
  grid[]
}

# ---- Main routine: download + aggregate one version -------------------------
process_version <- function(version, cfg) {
  out_dir    <- file.path(stations_dir, version)
  merged_csv <- file.path(stations_dir, sprintf("chirps_stations_%s.csv", version))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  grid <- build_grid(cfg, out_dir)
  
  # Download loop
  n <- nrow(grid)
  for (k in seq_len(n)) {
    r <- grid[k]
    if (file.exists(r$dest) && file.size(r$dest) > 0 && !overwrite) {
      grid[k, status := "skipped (exists)"]
    } else {
      ok <- download_one(r$url, r$dest, max_tries)
      grid[k, status := if (ok) "downloaded"
           else if (r$type %in% cfg$optional) "not available"
           else "failed"]
    }
    message(sprintf("[%s %4d/%d] %s -> %s", version, k, n, r$fname, grid$status[k]))
  }
  
  print(grid[, .N, by = .(type, status)])
  
  failed <- grid[status == "failed"]
  if (nrow(failed) > 0) {
    failed_csv <- file.path(stations_dir, sprintf("chirps_stations_%s_failed.csv", version))
    fwrite(failed[, .(year, month, url)], failed_csv)
    warning(sprintf("%s: %d file(s) failed; see %s", version, nrow(failed), failed_csv))
  }
  
  # Aggregate: read each file, filter to the target country, stack
  available <- grid[file.exists(dest)]
  merged <- rbindlist(
    lapply(seq_len(nrow(available)), function(k) {
      r  <- available[k]
      dt <- tryCatch(fread(r$dest, showProgress = FALSE), error = function(e) NULL)
      if (is.null(dt) || nrow(dt) == 0) return(NULL)
      dt <- filter_country(dt, target_country)
      if (is.null(dt) || nrow(dt) == 0) return(NULL)
      dt[, `:=`(year = r$year, month = r$month, file_type = r$type)]
      setcolorder(dt, c("year", "month", "file_type"))
      dt
    }),
    use.names = TRUE, fill = TRUE
  )
  
  if (nrow(merged) == 0) {
    warning(sprintf("%s: no rows for %s; nothing written.", version, target_country))
    return(invisible(merged))
  }
  
  setorder(merged, year, month)
  fwrite(merged, merged_csv)
  message(sprintf("%s: %d %s rows from %d files written to %s",
                  version, nrow(merged), target_country, nrow(available), merged_csv))
  invisible(merged)
}

# ---- Run: CHIRPS v3, then CHIRPS v2 -----------------------------------------
stations_v3 <- process_version("v3", versions$v3)
stations_v2 <- process_version("v2", versions$v2)