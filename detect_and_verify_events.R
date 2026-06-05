# ─────────────────────────────────────────────────────────────────────────────
#  Sargassum event detection and verification against on-site observations
#
#  Reads the *_coastal_positive_days.csv files (already filtered to slices
#  with >2 connected positive pixels and 5-day persistence), groups them into
#  discrete "events" using a coverage percentile threshold, and verifies
#  events against the all_observations.csv ground-truth dataset.
# ─────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(ggplot2)

# ───── Tunable parameters ────────────────────────────────────────────────────
COVERAGE_PERCENTILE <- 0.75   # threshold = this quantile of positive-day coverage_frac, per island
GAP_TOLERANCE       <- 10     # days; bridges a single cloud-obscured S2 acquisition
MIN_EVENT_DURATION  <- 5      # days; ≥2 Sentinel-2 acquisitions
OBS_BUFFER_DAYS     <- 7      # days; match window for ground-truth observations

# ───── Event detection ───────────────────────────────────────────────────────
detect_events <- function(pos_df,
                          percentile   = COVERAGE_PERCENTILE,
                          gap_tol      = GAP_TOLERANCE,
                          min_duration = MIN_EVENT_DURATION) {
  pos_df$date <- as.Date(pos_df$date)
  pos_df      <- pos_df[order(pos_df$date), ]

  # Coverage threshold from the distribution of positive-day coverage_frac
  threshold <- quantile(pos_df$coverage_frac, percentile, na.rm = TRUE)
  high      <- pos_df[pos_df$coverage_frac >= threshold, ]
  if (nrow(high) == 0) return(NULL)

  # Group consecutive high days into events, bridging gaps ≤ gap_tol
  gaps <- c(0, as.numeric(diff(high$date)))
  high$event_id <- cumsum(gaps > gap_tol) + 1L

  events <- high %>%
    group_by(event_id) %>%
    summarise(
      start          = min(date),
      end            = max(date),
      duration_days  = as.integer(max(date) - min(date)) + 1L,
      n_high_obs     = dplyr::n(),
      peak_coverage  = max(coverage_frac),
      mean_coverage  = mean(coverage_frac),
      peak_fai       = max(max_fai),
      .groups        = "drop"
    ) %>%
    filter(duration_days >= min_duration)

  attr(events, "threshold") <- threshold
  events
}

# ───── Observation matching ──────────────────────────────────────────────────
match_observations <- function(events, obs_df, buffer_days = OBS_BUFFER_DAYS) {
  # Strip tibble class and re-assert types defensively
  events       <- as.data.frame(events)
  events$start <- as.Date(events$start)
  events$end   <- as.Date(events$end)
  obs_df       <- as.data.frame(obs_df)
  obs_df$date  <- as.Date(obs_df$date)
  
  n <- nrow(events)
  n_obs     <- integer(n)
  first_obs <- as.Date(rep(NA_real_, n), origin = "1970-01-01")
  lag_days  <- rep(NA_integer_, n)
  
  valid <- !is.na(obs_df$date)
  for (i in seq_len(n)) {
    in_win <- valid &
      obs_df$date >= (events$start[i] - buffer_days) &
      obs_df$date <= (events$end[i]   + buffer_days)
    matched_dates <- obs_df$date[in_win]
    n_obs[i] <- length(matched_dates)
    if (length(matched_dates) > 0) {
      fo            <- min(matched_dates)
      first_obs[i]  <- fo
      lag_days[i]   <- as.integer(fo - events$start[i])
    }
  }
  
  events$n_observations <- n_obs
  events$first_obs_date <- first_obs
  events$obs_lag_days   <- lag_days
  events$verified       <- n_obs > 0
  events
}

# ───── Per-island wrapper ────────────────────────────────────────────────────
run_island <- function(island_name, pos_path, obs_all) {
  pos    <- read.csv(pos_path)
  events <- detect_events(pos)

  if (is.null(events) || nrow(events) == 0) {
    cat(sprintf("\n=== %s: no events detected ===\n", island_name))
    return(NULL)
  }

  obs_island      <- obs_all[obs_all$island == tolower(island_name), ]
  events          <- match_observations(events, obs_island)
  events$island   <- island_name

  thr <- attr(events, "threshold")
  cat(sprintf("\n=== %s ===\n", island_name))
  cat(sprintf("  positive days available:        %d\n", nrow(pos)))
  cat(sprintf("  threshold (P%.0f of coverage_frac): %.4f  (%.2f%%)\n",
              COVERAGE_PERCENTILE * 100, thr, thr * 100))
  cat(sprintf("  events detected:                %d\n", nrow(events)))
  cat(sprintf("  events verified by ≥1 obs:      %d (%.0f%%)\n",
              sum(events$verified), 100 * mean(events$verified)))
  cat(sprintf("  ground-truth observations:      %d\n", nrow(obs_island)))
  if (any(events$verified)) {
    cat(sprintf("  median obs lag from event start: %f d (negative = obs before satellite)\n",
                median(events$obs_lag_days, na.rm = TRUE)))
  }
  events
}

# ───── Reverse check: precision-side stat ────────────────────────────────────
reverse_check <- function(island_name, events, obs_all,
                          buffer_days = OBS_BUFFER_DAYS) {
  obs_island      <- as.data.frame(obs_all[obs_all$island == tolower(island_name), ])
  obs_island$date <- as.Date(obs_island$date)
  total_obs       <- nrow(obs_island)
  obs_island      <- obs_island[!is.na(obs_island$date), ]
  n_obs           <- nrow(obs_island)
  
  if (n_obs == 0) {
    cat(sprintf("  %s: no observations with valid dates\n", island_name))
    return(invisible())
  }
  if (is.null(events) || nrow(events) == 0) {
    cat(sprintf("  %s: no events detected; 0/%d observations explained\n",
                island_name, n_obs))
    return(invisible())
  }
  
  events       <- as.data.frame(events)
  events$start <- as.Date(events$start)
  events$end   <- as.Date(events$end)
  
  in_event <- logical(n_obs)
  for (i in seq_len(nrow(events))) {
    in_event <- in_event |
      (obs_island$date >= (events$start[i] - buffer_days) &
         obs_island$date <= (events$end[i]   + buffer_days))
  }
  
  cat(sprintf("  %s: %d / %d observations (%.0f%%) within ±%d d of a detected event  [%d obs dropped for unparseable dates]\n",
              island_name, sum(in_event), n_obs,
              100 * mean(in_event), buffer_days,
              total_obs - n_obs))
}
# ───── Visual diagnostic ─────────────────────────────────────────────────────
plot_island <- function(island_name, all_csv_path, events, obs_all) {
  # all_csv_path points to *_coastal_all.csv (every slice, not just positives) –
  # gives a complete timeline; if you don't have it, swap in positive_days CSV.
  all_df       <- read.csv(all_csv_path)
  all_df$date  <- as.Date(all_df$date)
  obs_island   <- obs_all[obs_all$island == tolower(island_name), ]

  p <- ggplot() +
    geom_col(data = all_df, aes(x = date, y = coverage_frac),
             fill = "#91D1C2", linewidth = 0.3)

  if (!is.null(events) && nrow(events) > 0) {
    p <- p + geom_rect(
      data = events,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = verified),
      alpha = 0.25, inherit.aes = FALSE
    ) + scale_fill_manual(
      values = c(`TRUE` = "#00A087", `FALSE` = "#DC0000"),
      name   = "verified by obs"
    )
  }

  p <- p +
    geom_point(data = obs_island, aes(x = date, y = 0),
               color = "black", size = 0.7, shape = 124) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title    = sprintf("%s: detected events vs ground-truth observations",
                         island_name),
      subtitle = sprintf("P%.0f threshold | min %d d | obs buffer ±%d d",
                         COVERAGE_PERCENTILE * 100,
                         MIN_EVENT_DURATION, OBS_BUFFER_DAYS),
      x        = NULL,
      y        = "coverage_frac"
    ) +
    theme_minimal()
  print(p)
}

# ───── Run ───────────────────────────────────────────────────────────────────
obs_all      <- read.csv("OBSdata/all_observations.csv")
obs_all$date <- as.Date(obs_all$date)

bonaire_events  <- run_island("Bonaire",
                              "sargassum_csvs/bonaire_coastal_positive_days.csv",  obs_all)
barbados_events <- run_island("Barbados",
                              "sargassum_csvs/barbados_coastal_positive_days.csv", obs_all)
aruba_events    <- run_island("Aruba",
                              "sargassum_csvs/aruba_coastal_positive_days.csv",    obs_all)

cat("\n=== Reverse check (observation → event) ===\n")
reverse_check("Bonaire",  bonaire_events,  obs_all)
reverse_check("Barbados", barbados_events, obs_all)
reverse_check("Aruba",    aruba_events,    obs_all)

events_all <- bind_rows(bonaire_events, barbados_events, aruba_events)
write.csv(events_all, "detected_events.csv", row.names = FALSE)

plot_island("Bonaire",  "sargassum_csvs/bonaire_coastal_all.csv",  bonaire_events,  obs_all)
plot_island("Barbados", "sargassum_csvs/barbados_coastal_all.csv", barbados_events, obs_all)
plot_island("Aruba",    "sargassum_csvs/aruba_coastal_all.csv",    aruba_events,    obs_all)
