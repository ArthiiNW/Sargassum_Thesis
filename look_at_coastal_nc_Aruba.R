library(ncdf4)
library(terra)
library(sf)
library(geodata)

# Let's also look at Aruba, which is known to have minimum
# Sargassum influx


nc_path     <- "data/Aruba.nc" 
varname_fai <- "nfai_max_isolated"
varname_obs <- "nfai_max_missing"
nc <- nc_open(nc_path)

lon <- ncvar_get(nc, "longitude")
lat <- ncvar_get(nc, "latitude")
nc_close(nc)

cat("lon: ", length(lon), " from ", min(lon), " to ", max(lon), "\n")
cat("lat: ", length(lat), " from ", min(lat), " to ", max(lat), "\n")
# Sanity check passed

# .　 . • ☆ . ° .• °:. *₊ ° . ☆MAP.　 . • ☆ . ° .• °:. *₊ ° . ☆

# --- Get Bonaire polygon from GADM ---
Aruba <- gadm(country = "ABW", level = 0, path = tempdir())
plot(Aruba)  # quick visual check

# --- Build 1 km seaward buffer in a metric CRS (UTM 19N) ---
utm <- "EPSG:32619"
Aruba_utm <- project(Aruba, utm)
buffer_utm  <- buffer(Aruba_utm, width = 1000)   # 1 km outward
seaward_utm <- erase(buffer_utm, Aruba_utm)      # ring only, drops land
seaward     <- project(seaward_utm, "EPSG:4326")

# --- Determine which grid cells fall inside the buffer ---
# Build coordinates in the same order R fills a matrix(values, n_lon, n_lat):
# columns vary slowest, so lat varies slowest in `each =`.
coords <- cbind(
  lon = rep(lon, times = length(lat)),
  lat = rep(lat, each  = length(lon))
)
pts <- vect(coords, type = "points", crs = "EPSG:4326")

hits <- relate(pts, seaward, "intersects")[, 1]
buffer_mask <- matrix(hits, nrow = length(lon), ncol = length(lat))

cat("Cells inside 1 km buffer:", sum(buffer_mask), "\n")
# 1477 for Bonaire vs 1318 for Barbados
# Sanity-check the mask visually (Bonaire-shaped ring should appear)
image(lon, lat, buffer_mask, asp = 1,
      main = "1 km seaward buffer (TRUE cells)")
plot(Aruba, add = TRUE, border = "red")
# sanity-check passed

# .　 . • ☆ . ° .• °:. *₊ ° . ☆NC.　 . • ☆ . ° .• °:. *₊ ° . ☆

nc <- nc_open(nc_path)
# Find a slice with FAI activity
n_time <- nc$dim$time$len

n_buffer_cells <- sum(buffer_mask)

# Decode time
time_raw   <- ncvar_get(nc, "time")
time_units <- ncatt_get(nc, "time", "units")$value
parse_time <- function(tr, tu) {
  parts <- strsplit(tu, " since ")[[1]]
  mult  <- switch(tolower(parts[1]), days=86400, hours=3600, minutes=60, seconds=1)
  as.POSIXct(parts[2], tz="UTC") + tr * mult
}
dates <- as.Date(parse_time(time_raw, time_units))

max_fai      <- numeric(n_time)
n_positive   <- integer(n_time)
n_observed   <- integer(n_time)   # observed within buffer

for (t in seq_len(n_time)) {
  fai_slice <- ncvar_get(nc, varname_fai, start=c(1,1,t), count=c(-1,-1,1))
  
  # Primary validity mask: FAI is not NA. Restrict to coastal buffer.
  in_scope <- buffer_mask & !is.na(fai_slice)
  
  fai_in_scope <- fai_slice
  fai_in_scope[!in_scope] <- NA
  
  max_fai[t]    <- suppressWarnings(max(fai_in_scope, na.rm = TRUE))
  n_positive[t] <- sum(fai_in_scope > 0, na.rm = TRUE)
  n_observed[t] <- sum(in_scope)
  
  if (t %% 200 == 0) cat("processed", t, "/", n_time, "\n")
}
max_fai[is.infinite(max_fai)] <- NA
nc_close(nc)

results <- data.frame(
  slice_index   = seq_len(n_time),
  date          = dates,
  max_fai       = max_fai,
  n_positive    = n_positive,
  n_observed    = n_observed,
  coverage_frac = n_observed / n_buffer_cells # how much of the coastal area are covered?
)

positive_days <- subset(
  results,
  !is.na(max_fai) & max_fai > 0 # does not make sense to have a coverage filter: very small coverage
)

cat("Coastal positive-FAI days:", nrow(positive_days), "\n")
cat("Total days with any coastal observation:", sum(results$n_observed > 0), "\n")
cat("Median coverage on observed days:",
    median(results$coverage_frac[results$n_observed > 0]), "\n")

# write.csv(positive_days, "bonaire_coastal_positive_days.csv", row.names=FALSE)
# write.csv(results,       "bonaire_coastal_summary.csv",       row.names=FALSE)

# .　 . • ☆ . ° .• °:. *₊ ° . ☆Look at data.　 . • ☆ . ° .• °:. *₊ ° . ☆

library(ggplot2)

positive_days$year       <- as.integer(format(positive_days$date, "%Y"))
positive_days$month      <- as.integer(format(positive_days$date, "%m"))
positive_days$month_name <- factor(format(positive_days$date, "%b"),
                                   levels = month.abb)


# (a) Timeline: each positive day as a tick, colored by intensity
p_timeline <- ggplot(results, aes(x = date, y = coverage_frac)) +
  geom_area(fill = "darkorange", alpha = 0.4) +
  geom_line(color = "darkorange", linewidth = 0.3) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Sargassum coverage in Aruba 5 km coastal buffer",
       x = NULL, y = "fraction of buffer cells with positive FAI") +
  theme_minimal()
print(p_timeline)

# Same pattern with Bonaire

results$smoothed <- zoo::rollmean(results$coverage_frac, k = 14,
                                  fill = NA, align = "center")

p_timeline + geom_line(aes(y = results$smoothed), color = "firebrick", linewidth = 0.7)

results$year <- as.integer(format(results$date, "%Y"))
results$doy  <- as.integer(format(results$date, "%j"))

p_calendar <- ggplot(results, aes(x = doy, y = factor(year),
                                  fill = coverage_frac)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", direction = -1, na.value = "grey95",
                       labels = scales::percent,
                       name = "sargassum\ncoverage") +
  scale_x_continuous(
    breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
    labels = month.abb, expand = c(0, 0)) +
  labs(title = "Daily sargassum coverage in Bonaire coastal buffer",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank())
print(p_calendar)

# (b) Year-month heatmap: shows seasonality + interannual variation at a glance
ym_counts <- as.data.frame(table(year = positive_days$year,
                                 month = positive_days$month_name))
p_heatmap <- ggplot(ym_counts, aes(x = month, y = factor(year), fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "mako", direction = -1,
                       name = "positive days") +
  labs(title = "Positive-FAI days by year and month",
       x = NULL, y = NULL) +
  theme_minimal()

# (c) Bar chart of positive days per year
p_year <- ggplot(positive_days, aes(x = factor(year))) +
  geom_bar(fill = "steelblue") +
  labs(title = "Positive-FAI days per year", x = "year", y = "days") +
  theme_minimal()

print(p_timeline); print(p_heatmap); print(p_year)







library(ncdf4)
library(terra)

nc <- nc_open(nc_path)

# Accumulators (matching the lon × lat grid orientation)
freq_pos    <- matrix(0L, nrow = length(lon), ncol = length(lat))   # count of positive observations
freq_obs    <- matrix(0L, nrow = length(lon), ncol = length(lat))   # count of valid observations
sum_fai_pos <- matrix(0,  nrow = length(lon), ncol = length(lat))   # for mean positive intensity

for (i in seq_len(nrow(positive_days))) {
  t   <- positive_days$slice_index[i]
  fai <- ncvar_get(nc, varname_fai, start = c(1, 1, t), count = c(-1, -1, 1))
  
  observed <- !is.na(fai) & buffer_mask
  positive <- observed & fai > 0
  
  freq_obs    <- freq_obs    + observed
  freq_pos    <- freq_pos    + positive
  sum_fai_pos[positive] <- sum_fai_pos[positive] + fai[positive]
  
  if (i %% 50 == 0) cat("processed", i, "/", nrow(positive_days), "\n")
}
nc_close(nc)

# Frequency of positive detections (raw count and as a rate per observation)
rate_pos <- ifelse(freq_obs > 0, freq_pos / freq_obs, NA)
mean_intensity <- ifelse(freq_pos > 0, sum_fai_pos / freq_pos, NA)

# Mask to coastal buffer for cleaner plotting
freq_pos_masked       <- freq_pos
rate_pos_masked       <- rate_pos
mean_intensity_masked <- mean_intensity
freq_pos_masked[!buffer_mask | freq_pos == 0]       <- NA
rate_pos_masked[!buffer_mask]                       <- NA
mean_intensity_masked[!buffer_mask]                 <- NA

# Build SpatRasters for nice plotting with terra
make_rast <- function(mat) {
  # terra expects rows = lat (top-to-bottom), cols = lon (left-to-right).
  # mat is [lon, lat] with lat ascending → transpose, then flip vertically.
  r <- rast(t(mat)[length(lat):1, ],
            extent = ext(min(lon), max(lon), min(lat), max(lat)),
            crs = "EPSG:4326")
  r
}

r_count     <- make_rast(freq_pos_masked)
r_rate      <- make_rast(rate_pos_masked)
r_intensity <- make_rast(mean_intensity_masked)

# Plot
par(mfrow = c(1, 3), mar = c(3, 3, 3, 5))

plot(r_count, main = "Positive detections (count)",
     col = hcl.colors(20, "YlOrRd"))
plot(Aruba, add = TRUE, border = "black", lwd = 1.2)

plot(r_rate, main = "Positive rate (positives / observations)",
     col = hcl.colors(20, "YlOrRd"))
plot(Aruba, add = TRUE, border = "black", lwd = 1.2)

plot(r_intensity, main = "Mean FAI when positive",
     col = hcl.colors(20, "YlOrRd"))
plot(Aruba, add = TRUE, border = "black", lwd = 1.2)

par(mfrow = c(1, 1))

