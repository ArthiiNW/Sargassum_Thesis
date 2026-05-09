library(ncdf4)

nc_path <- "data/Isolate_max_noob_FAI.nc" 
varname <- "nfai_max_isolated"

nc <- nc_open(nc_path)
n_time <- nc$dim$time$len

# --- Decode time axis ---
time_raw   <- ncvar_get(nc, "time")
time_units <- ncatt_get(nc, "time", "units")$value
cat("time units:", time_units, "\n")   # sanity check

parse_time <- function(time_raw, time_units) {
  parts  <- strsplit(time_units, " since ")[[1]]
  unit   <- tolower(parts[1])
  origin <- as.POSIXct(parts[2], tz = "UTC")
  mult   <- switch(unit,
                   "days"    = 86400,
                   "hours"   = 3600,
                   "minutes" = 60,
                   "seconds" = 1,
                   stop("Unknown time unit: ", unit))
  origin + time_raw * mult
}
dates <- as.Date(parse_time(time_raw, time_units))

# --- Chunked pass: one time slice at a time ---
max_fai    <- numeric(n_time)
n_positive <- integer(n_time)
n_valid    <- integer(n_time)

for (t in seq_len(n_time)) {
  slice <- ncvar_get(nc, varname,
                     start = c(1, 1, t),
                     count = c(-1, -1, 1))
  max_fai[t]    <- suppressWarnings(max(slice, na.rm = TRUE))
  n_positive[t] <- sum(slice > 0, na.rm = TRUE)
  n_valid[t]    <- sum(!is.na(slice))
  if (t %% 200 == 0) cat("processed", t, "/", n_time, "\n")
}

# All-NaN slices produce -Inf from max(); flag as NA
max_fai[is.infinite(max_fai)] <- NA

nc_close(nc)

# --- Build full results table (keep slice_index so you can map back) ---
results <- data.frame(
  slice_index = seq_len(n_time),
  date        = dates,
  max_fai     = max_fai,
  n_positive  = n_positive,
  n_valid     = n_valid
)

# --- Filter to positive-FAI days ---
positive_days <- subset(results, !is.na(max_fai) & max_fai > 0)

cat("Positive-FAI days found:", nrow(positive_days),
    "out of", sum(!is.na(max_fai)), "valid days\n")

get_slice <- function(nc_path, slice_index, varname = "nfai_max_isolated") {
  nc <- nc_open(nc_path)
  on.exit(nc_close(nc))
  ncvar_get(nc, varname,
            start = c(1, 1, slice_index),
            count = c(-1, -1, 1))
}

# Example: pull the map for the first positive-FAI day
i      <- positive_days$slice_index[7]
d      <- positive_days$date[7]
map_ij <- get_slice(nc_path, i)

image(map_ij, main = paste("NFAI on", d), col = hcl.colors(20, "YlOrRd"))

# Or pull a map by date:
target_date <- as.Date("2018-08-15")
i <- results$slice_index[results$date == target_date]
map_ij <- get_slice(nc_path, i)