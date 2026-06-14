# =============================================================
# Diagnostic: observation effort by month-of-year
# -------------------------------------------------------------
# Purpose: assess whether clear-look satellite acquisitions over the
#   coastal buffers are seasonally structured. If wet-season months are
#   systematically under-observed, then recording an unobserved month as
#   zero sargassum (missing-as-false) injects a seasonally structured
#   (MNAR) misclassification bias rather than mean-zero noise -- the very
#   problem flagged in the misclassification sub-section.
#
# Input : <island>_coastal_all.csv  (daily grid; one row per calendar day.
#   n_observed = count of cloud-free pixels in the coastal buffer that day;
#   n_observed == 0 means no clear look, i.e. no pass or fully obscured.)
#
# Output: fig_observation_effort_by_month.png  (primary: acquisition rate)
#         fig_view_completeness_by_month.png   (optional: masking signal)
#         + a tidy summary table printed to the console.
# =============================================================


# ---- 0. Packages -------------------------------------------------------
library(dplyr)
library(ggplot2)


# ---- 1. Read the three daily coastal series and stack them -------------
data_dir <- "sargassum_csvs"

read_island <- function(file, island_name) {
  df <- read.csv(file.path(data_dir, file), stringsAsFactors = FALSE)
  df$island <- island_name
  df
}

coastal <- bind_rows(
  read_island("bonaire_coastal_all.csv",  "Bonaire"),
  read_island("aruba_coastal_all.csv",    "Aruba"),
  read_island("barbados_coastal_all.csv", "Barbados")
)

# Calendar fields from the date string.
coastal$date  <- as.Date(coastal$date)
coastal$month <- as.integer(format(coastal$date, "%m"))

# A "clear look" = at least one cloud-free pixel in the buffer that day.
coastal$observed <- coastal$n_observed > 0


# ---- 2. Aggregate by island x month-of-year ----------------------------
# NOTE on fair comparison across months: the series does not cover an equal
# number of years per calendar month (Dec begins mid-Dec 2020; Jan-Mar
# include 2026; Apr-Nov stop in 2025). A raw count of clear-look days would
# therefore be biased by how many years each month happens to be sampled.
# We instead use the ACQUISITION RATE = clear-look days / calendar days in
# the grid for that month-of-year, which is invariant to unequal years.
effort <- coastal |>
  group_by(island, month) |>
  summarise(
    calendar_days = n(),                        # daily rows available
    acq_days      = sum(observed),              # days with a clear look
    acq_rate      = acq_days / calendar_days,   # detection opportunity
    mean_nobs     = mean(n_observed[observed]), # view completeness | look
    .groups = "drop"
  )

# Plotting order.
effort$month_lab <- factor(month.abb[effort$month], levels = month.abb)
effort$island    <- factor(effort$island,
                           levels = c("Bonaire", "Aruba", "Barbados"))


# ---- 3. Primary diagnostic: acquisition rate by month-of-year ----------
# Caribbean wet season (~Jun-Nov) shaded for context. Delete the geom_rect
# block if you prefer an unshaded figure. The shaded band is numeric on a
# discrete x-axis: positions Jun=6 .. Nov=11, so 5.5-11.5 covers them.
wet_band <- data.frame(xmin = 5.5, xmax = 11.5, ymin = -Inf, ymax = Inf)

p_rate <- ggplot(effort, aes(x = month_lab, y = acq_rate)) +
  geom_rect(data = wet_band, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey85", alpha = 0.5) +
  geom_col(fill = "#2c7fb8", width = 0.7) +
  facet_wrap(~ island, ncol = 1, scales = "free_y") +
  labs(
    title    = "Observation effort by month-of-year",
    subtitle = paste0("Fraction of calendar days yielding a clear coastal ",
                      "look (n_observed > 0).\n",
                      "Shaded band = approximate Caribbean wet season (Jun-Nov)."),
    x = NULL,
    y = "Acquisition rate (clear-look days / calendar days)"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("fig_observation_effort_by_month.png", p_rate,
       width = 7, height = 7, dpi = 300)


# ---- 4. (Optional) View completeness: mean cloud-free pixels per look ---
# Speaks to the INTENSITY-suppression mechanism rather than the frequency
# one: if a clear look in the wet season still reveals fewer usable pixels,
# coverage/intensity is mechanically depressed in exactly the months that
# overlap the sargassum season. Note the cross-island level differences
# reflect buffer size / UTM zone, so compare the seasonal *shape* within an
# island, not the absolute heights between islands (hence free y-scales).
p_complete <- ggplot(effort, aes(x = month_lab, y = mean_nobs, group = island)) +
  geom_rect(data = wet_band, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey85", alpha = 0.5) +
  geom_line(colour = "#d95f0e", linewidth = 0.8) +
  geom_point(colour = "#d95f0e", size = 1.8) +
  facet_wrap(~ island, ncol = 1, scales = "free_y") +
  labs(
    title    = "View completeness by month-of-year",
    subtitle = "Mean cloud-free pixels per clear-look day (conditional on a look).",
    x = NULL,
    y = "Mean n_observed | clear look"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("fig_view_completeness_by_month.png", p_complete,
       width = 7, height = 7, dpi = 300)


# ---- 5. Print the underlying table (handy for the thesis text) ---------
print(
  effort[order(effort$island, effort$month),
         c("island", "month", "calendar_days", "acq_days", "acq_rate", "mean_nobs")],
  row.names = FALSE, digits = 3
)
