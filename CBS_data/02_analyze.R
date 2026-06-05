# ============================================================================
# 02_analyze.R
#
# Linear analysis: sargassum × tourism on Bonaire.
# Run section by section in RStudio (Ctrl+Enter or Cmd+Enter).
# Every intermediate object stays in the global environment so you can
# View() it, head() it, summarise it, plot it -- whatever you want.
#
# Prerequisites:
#   - bonaire_coastal_positive_days.csv in the working directory
#   - cbs_raw/<id>.csv from running 01_fetch_cbs.R
#
# Sargassum variables used (only those whose interpretation is supported
# by the FAI literature; max_fai is excluded as it's a relative contrast
# measure without absolute meaning):
#   any_sargassum         -- was sargassum detected this month?
#   positive_days         -- # days with FAI-positive coastal pixels (severity)
#   total_positive_pixels -- summed positive pixels (landing-volume proxy)
# ============================================================================

# ---- 0. Setup -------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(ggplot2)
library(sandwich)
library(lmtest)
library(broom)

dir.create("output",         showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 11))

# We'll accumulate one-row summaries from each regression into this list,
# then bind them at the very end into a single results table.
results <- list()


# ============================================================================
# 1. SARGASSUM: raw FAI observations -> monthly severity series
# ============================================================================

# Load raw data
sarg_raw <- read_csv("sargassum_csvs/bonaire_coastal_positive_days.csv", show_col_types = FALSE)
glimpse(sarg_raw)
# NOTE: this CSV contains ONLY days where >=1 coastal pixel was FAI-positive.
# Months with no detection are absent from the file and must be inserted
# as zeros below, otherwise they would silently drop out of regressions.

sarg_raw <- sarg_raw |> mutate(date = as.Date(date))
range(sarg_raw$date)

# Build complete monthly index spanning the observation period
months_idx <- tibble(
  month_start = seq(floor_date(min(sarg_raw$date), "month"),
                    floor_date(max(sarg_raw$date), "month"),
                    by = "month")
) |>
  mutate(year = year(month_start), month = month(month_start))

head(months_idx); tail(months_idx)

# Aggregate detected-sargassum days per (year, month)
sarg_agg <- sarg_raw |>
  group_by(year, month) |>
  summarise(
    positive_days         = n_distinct(date),
    total_positive_pixels = sum(n_positive),
    .groups = "drop"
  )

# Merge into the complete index; missing months -> 0
sarg_monthly <- months_idx |>
  left_join(sarg_agg, by = c("year", "month")) |>
  mutate(
    positive_days         = replace_na(positive_days, 0L),
    total_positive_pixels = replace_na(total_positive_pixels, 0L),
    any_sargassum         = as.integer(positive_days > 0L)
  )

# Inspect
head(sarg_monthly, 14)
sarg_monthly |> count(any_sargassum)
sarg_monthly |> group_by(year) |>
  summarise(months_with_sarg = sum(any_sargassum),
            total_pos_days   = sum(positive_days),
            total_pos_pixels = sum(total_positive_pixels))

# Quick visual sanity check
ggplot(sarg_monthly, aes(month_start, positive_days)) +
  geom_col(fill = "#E48F4A") +
  labs(title = "Sargassum positive days per month — Bonaire coast",
       x = NULL, y = "positive days / month")

write_csv(sarg_monthly, "output/sargassum_monthly.csv")

# Also build a yearly roll-up for the value_per_sector dataset
sarg_yearly <- sarg_monthly |>
  group_by(year) |>
  summarise(
    positive_days         = sum(positive_days),
    positive_months       = sum(any_sargassum),
    total_positive_pixels = sum(total_positive_pixels)
  ) |>
  mutate(any_sargassum = as.integer(positive_days > 0L))
sarg_yearly


# ============================================================================
# 2. TOURISM BY PLANE  (83104NED) -- monthly, with region dimension
# ============================================================================

plane_raw <- read_csv("cbs_raw/83104NED.csv", show_col_types = FALSE)
glimpse(plane_raw)
# Note the column structure: Perioden_Date and Perioden_freq are added by
# cbs_add_date_column(); RegioS_label by cbs_add_label_columns().

# Filter to Bonaire (using the human-readable label) and extract year/month
plane_bo <- plane_raw |>
  filter(RegioS_label == "Bonaire") |>
  mutate(
    year         = year(as.Date(Perioden_Date)),
    month        = month(as.Date(Perioden_Date)),
    month_start  = floor_date(as.Date(Perioden_Date), "month")
  ) |>
  select(year, month, month_start, Perioden_freq, Aantal_1, Gemiddelde_2)
head(plane_bo)
# Aantal_1     = number of tourists by plane
# Gemiddelde_2 = average something (check meta file for unit)

# Join with sargassum
plane_m <- inner_join(sarg_monthly, plane_bo, by = c("year", "month", "month_start"))
glimpse(plane_m)

# Plot
ggplot(plane_m, aes(month_start)) +
  geom_col(aes(y = positive_days * max(Aantal_1, na.rm = TRUE) /
                   max(positive_days, 1)), fill = "#E48F4A", alpha = 0.4) +
  geom_line(aes(y = Aantal_1), colour = "#2E5C8A", linewidth = 0.9) +
  labs(title = "Tourism by plane vs sargassum (83104NED)",
       y = "Aantal_1 (tourists by plane)", x = NULL) -> p_plane
print(p_plane)
ggsave("output/figures/plane_timeseries.png", p_plane, width = 9, height = 4, dpi = 120)

# Regression: tourists ~ positive_days + month FE + year trend, HAC SE
plane_m <- plane_m |>
  mutate(month_fac = factor(month), year_trend = year - min(year))
m_plane <- lm(Aantal_1 ~ positive_days + month_fac + year_trend, data = plane_m)
summary(m_plane)            # classical SEs
vc_plane <- NeweyWest(m_plane, lag = 3, prewhite = FALSE)
coeftest(m_plane, vcov = vc_plane)  # HAC SEs (the ones to report)

# Tidy result row for the summary table
res_plane <- tidy(coeftest(m_plane, vcov = vc_plane), conf.int = TRUE,
                  conf.level = 0.95) |>
  filter(term == "positive_days") |>
  mutate(dataset = "tourism_by_plane", series = "Total", measure = "Aantal_1",
         n_obs = nrow(plane_m), freq = "M")
res_plane
results[["plane_Aantal"]] <- res_plane


# ============================================================================
# 3. TOURISM BY NATIONALITY  (83191NED) -- monthly, nationality sub-dimension
# ============================================================================

nat_raw <- read_csv("cbs_raw/83191NED.csv", show_col_types = FALSE)
glimpse(nat_raw)

# All rows are already Bonaire-only for this table, but be safe
nat_bo <- nat_raw |>
  filter(RegioS_label == "Bonaire") |>
  mutate(year = year(as.Date(Perioden_Date)),
         month = month(as.Date(Perioden_Date)),
         month_start = floor_date(as.Date(Perioden_Date), "month")) |>
  select(year, month, month_start, Nationaliteit_label, Aantal_1)
head(nat_bo)
unique(nat_bo$Nationaliteit_label)

# Total across nationalities
nat_total <- nat_bo |>
  group_by(year, month, month_start) |>
  summarise(Aantal_1 = sum(Aantal_1, na.rm = TRUE), .groups = "drop")
head(nat_total)

# Wide form by nationality (one column per nationality) for plotting
nat_wide <- nat_bo |>
  pivot_wider(names_from = Nationaliteit_label, values_from = Aantal_1,
              values_fn = sum) |>
  arrange(month_start)
head(nat_wide)

# Merge total with sargassum
nat_m <- inner_join(sarg_monthly, nat_total, by = c("year", "month", "month_start"))

# Plot all nationalities + total
nat_long <- bind_rows(
  nat_bo |> rename(value = Aantal_1, group = Nationaliteit_label),
  nat_total |> mutate(group = "Total") |> rename(value = Aantal_1)
)
ggplot(nat_long, aes(month_start, value, colour = group)) +
  geom_line() +
  labs(title = "Tourism by nationality (83191NED) — Bonaire",
       y = "monthly arrivals", x = NULL, colour = NULL) -> p_nat
print(p_nat)
ggsave("output/figures/nationality_timeseries.png", p_nat, width = 9, height = 4, dpi = 120)

# Regression on TOTAL
nat_m <- nat_m |> mutate(month_fac = factor(month), year_trend = year - min(year))
m_nat_total <- lm(Aantal_1 ~ positive_days + month_fac + year_trend, data = nat_m)
vc_nat <- NeweyWest(m_nat_total, lag = 3, prewhite = FALSE)
coeftest(m_nat_total, vcov = vc_nat)

results[["nationality_Total"]] <- tidy(coeftest(m_nat_total, vcov = vc_nat),
                                        conf.int = TRUE, conf.level = 0.95) |>
  filter(term == "positive_days") |>
  mutate(dataset = "tourism_by_nationality", series = "Total",
         measure = "Aantal_1", n_obs = nrow(nat_m), freq = "M")

# Also regress per-nationality so you can see if e.g. US vs NL tourists react differently
for (nat in unique(nat_bo$Nationaliteit_label)) {
  d_nat <- nat_bo |> filter(Nationaliteit_label == nat) |>
    inner_join(sarg_monthly, by = c("year", "month", "month_start")) |>
    mutate(month_fac = factor(month), year_trend = year - min(year))
  if (nrow(d_nat) < 12) next
  m_n <- lm(Aantal_1 ~ positive_days + month_fac + year_trend, data = d_nat)
  vc_n <- NeweyWest(m_n, lag = 3, prewhite = FALSE)
  cat(sprintf("\n--- Nationality: %s (n=%d) ---\n", nat, nrow(d_nat)))
  print(coeftest(m_n, vcov = vc_n)["positive_days", , drop = FALSE])
  results[[paste0("nationality_", nat)]] <- tidy(coeftest(m_n, vcov = vc_n),
                                                  conf.int = TRUE) |>
    filter(term == "positive_days") |>
    mutate(dataset = "tourism_by_nationality", series = nat,
           measure = "Aantal_1", n_obs = nrow(d_nat), freq = "M")
}


# ============================================================================
# 4. YACHTS & MOORING DAYS  (85015NED) -- monthly, Bonaire-only
# ============================================================================

yacht_raw <- read_csv("cbs_raw/85015NED.csv", show_col_types = FALSE)
glimpse(yacht_raw)
# No RegioS column -- this table is Bonaire-only by construction

yacht_bo <- yacht_raw |>
  mutate(year = year(as.Date(Perioden_Date)),
         month = month(as.Date(Perioden_Date)),
         month_start = floor_date(as.Date(Perioden_Date), "month")) |>
  select(year, month, month_start, AantalJachten_1, Aanmeerdagen_2)
head(yacht_bo)

yacht_m <- inner_join(sarg_monthly, yacht_bo, by = c("year", "month", "month_start"))

ggplot(yacht_m, aes(month_start)) +
  geom_col(aes(y = positive_days * max(AantalJachten_1, na.rm = TRUE) /
                   max(positive_days, 1)), fill = "#E48F4A", alpha = 0.4) +
  geom_line(aes(y = AantalJachten_1), colour = "#2E5C8A", linewidth = 0.9) +
  labs(title = "Yachts vs sargassum (85015NED)",
       y = "yachts / month", x = NULL) -> p_yacht
print(p_yacht)
ggsave("output/figures/yachts_timeseries.png", p_yacht, width = 9, height = 4, dpi = 120)

yacht_m <- yacht_m |> mutate(month_fac = factor(month), year_trend = year - min(year))
m_yacht <- lm(AantalJachten_1 ~ positive_days + month_fac + year_trend, data = yacht_m)
vc_yacht <- NeweyWest(m_yacht, lag = 3, prewhite = FALSE)
coeftest(m_yacht, vcov = vc_yacht)
results[["yachts_AantalJachten"]] <- tidy(coeftest(m_yacht, vcov = vc_yacht),
                                          conf.int = TRUE) |>
  filter(term == "positive_days") |>
  mutate(dataset = "yachts_mooring", series = "Total",
         measure = "AantalJachten_1", n_obs = nrow(yacht_m), freq = "M")

# Same for mooring days
m_moor <- lm(Aanmeerdagen_2 ~ positive_days + month_fac + year_trend, data = yacht_m)
vc_moor <- NeweyWest(m_moor, lag = 3, prewhite = FALSE)
coeftest(m_moor, vcov = vc_moor)
results[["yachts_Aanmeerdagen"]] <- tidy(coeftest(m_moor, vcov = vc_moor),
                                          conf.int = TRUE) |>
  filter(term == "positive_days") |>
  mutate(dataset = "yachts_mooring", series = "Total",
         measure = "Aanmeerdagen_2", n_obs = nrow(yacht_m), freq = "M")


# ============================================================================
# 5. AVIATION  (82332NED) -- monthly, with region dimension
# ============================================================================

avi_raw <- read_csv("cbs_raw/82332NED.csv", show_col_types = FALSE)
glimpse(avi_raw)

avi_bo <- avi_raw |>
  filter(RegioS_label == "Bonaire") |>
  mutate(year = year(as.Date(Perioden_Date)),
         month = month(as.Date(Perioden_Date)),
         month_start = floor_date(as.Date(Perioden_Date), "month")) |>
  select(year, month, month_start, Passagiers_1, Vluchten_2)
head(avi_bo)

avi_m <- inner_join(sarg_monthly, avi_bo, by = c("year", "month", "month_start")) |>
  mutate(month_fac = factor(month), year_trend = year - min(year))

m_pax <- lm(Passagiers_1 ~ positive_days + month_fac + year_trend, data = avi_m)
vc_pax <- NeweyWest(m_pax, lag = 3, prewhite = FALSE)
coeftest(m_pax, vcov = vc_pax)
results[["aviation_Passagiers"]] <- tidy(coeftest(m_pax, vcov = vc_pax),
                                          conf.int = TRUE) |>
  filter(term == "positive_days") |>
  mutate(dataset = "aviation", series = "Total",
         measure = "Passagiers_1", n_obs = nrow(avi_m), freq = "M")

m_flt <- lm(Vluchten_2 ~ positive_days + month_fac + year_trend, data = avi_m)
vc_flt <- NeweyWest(m_flt, lag = 3, prewhite = FALSE)
coeftest(m_flt, vcov = vc_flt)
results[["aviation_Vluchten"]] <- tidy(coeftest(m_flt, vcov = vc_flt),
                                        conf.int = TRUE) |>
  filter(term == "positive_days") |>
  mutate(dataset = "aviation", series = "Total",
         measure = "Vluchten_2", n_obs = nrow(avi_m), freq = "M")


# ============================================================================
# 6. CRUISE PASSENGERS  (85007NED) -- monthly, Bonaire-only
# ============================================================================

cru_raw <- read_csv("cbs_raw/85007NED.csv", show_col_types = FALSE)
glimpse(cru_raw)

cru_bo <- cru_raw |>
  mutate(year = year(as.Date(Perioden_Date)),
         month = month(as.Date(Perioden_Date)),
         month_start = floor_date(as.Date(Perioden_Date), "month")) |>
  select(year, month, month_start, Aantal_1, AantalSchepen_2)
head(cru_bo)
# Aantal_1        = cruise passengers
# AantalSchepen_2 = number of cruise ships

cru_m <- inner_join(sarg_monthly, cru_bo, by = c("year", "month", "month_start")) |>
  mutate(month_fac = factor(month), year_trend = year - min(year))

ggplot(cru_m, aes(month_start)) +
  geom_col(aes(y = positive_days * max(Aantal_1, na.rm = TRUE) /
                   max(positive_days, 1)), fill = "#E48F4A", alpha = 0.4) +
  geom_line(aes(y = Aantal_1), colour = "#2E5C8A", linewidth = 0.9) +
  labs(title = "Cruise passengers vs sargassum (85007NED)",
       y = "passengers / month", x = NULL) -> p_cru
print(p_cru)
ggsave("output/figures/cruise_timeseries.png", p_cru, width = 9, height = 4, dpi = 120)

m_cru <- lm(Aantal_1 ~ positive_days + month_fac + year_trend, data = cru_m)
vc_cru <- NeweyWest(m_cru, lag = 3, prewhite = FALSE)
coeftest(m_cru, vcov = vc_cru)
results[["cruise_Aantal"]] <- tidy(coeftest(m_cru, vcov = vc_cru),
                                    conf.int = TRUE) |>
  filter(term == "positive_days") |>
  mutate(dataset = "cruise_passengers", series = "Total",
         measure = "Aantal_1", n_obs = nrow(cru_m), freq = "M")

m_ships <- lm(AantalSchepen_2 ~ positive_days + month_fac + year_trend, data = cru_m)
vc_ships <- NeweyWest(m_ships, lag = 3, prewhite = FALSE)
coeftest(m_ships, vcov = vc_ships)
results[["cruise_AantalSchepen"]] <- tidy(coeftest(m_ships, vcov = vc_ships),
                                           conf.int = TRUE) |>
  filter(term == "positive_days") |>
  mutate(dataset = "cruise_passengers", series = "Total",
         measure = "AantalSchepen_2", n_obs = nrow(cru_m), freq = "M")


# ============================================================================
# 7. VALUE PER SECTOR  (84769NED) -- YEARLY, with region + sector dimensions
# ============================================================================
# This table is annual, so we join on year only and lose all within-year
# variation. With ~7 post-2020 years of overlap the result is exploratory
# at best -- no month FE is possible.

sec_raw <- read_csv("cbs_raw/84769NED.csv", show_col_types = FALSE)
glimpse(sec_raw)

sec_bo <- sec_raw |>
  filter(RegioS_label == "Bonaire") |>
  mutate(year = year(as.Date(Perioden_Date))) |>
  select(year, Bedrijfstak_label, ToegevoegdeWaarde_1)
head(sec_bo)
unique(sec_bo$Bedrijfstak_label)

# Join with yearly sargassum
sec_y <- inner_join(sarg_yearly, sec_bo, by = "year")
sec_y |> arrange(Bedrijfstak_label, year)

ggplot(sec_y, aes(year, ToegevoegdeWaarde_1, colour = Bedrijfstak_label)) +
  geom_line() + geom_point() +
  labs(title = "Value-added by sector (84769NED) — Bonaire",
       y = "toegevoegde waarde", x = NULL, colour = "sector") -> p_sec
print(p_sec)
ggsave("output/figures/sector_timeseries.png", p_sec, width = 9, height = 4.5, dpi = 120)

# Per-sector regression on yearly positive_days (no month FE possible)
for (sect in unique(sec_y$Bedrijfstak_label)) {
  d_s <- sec_y |> filter(Bedrijfstak_label == sect) |>
    mutate(year_trend = year - min(year))
  if (nrow(d_s) < 5) next
  m_s <- lm(ToegevoegdeWaarde_1 ~ positive_days + year_trend, data = d_s)
  cat(sprintf("\n--- Sector: %s (n=%d) ---\n", sect, nrow(d_s)))
  print(summary(m_s)$coefficients["positive_days", , drop = FALSE])
  results[[paste0("sector_", sect)]] <- tidy(m_s, conf.int = TRUE) |>
    filter(term == "positive_days") |>
    mutate(dataset = "value_per_sector", series = sect,
           measure = "ToegevoegdeWaarde_1", n_obs = nrow(d_s), freq = "Y")
}


# ============================================================================
# 8. COMBINED RESULTS TABLE
# ============================================================================

summary_tbl <- bind_rows(results) |>
  select(dataset, series, measure, freq, n_obs,
         coef = estimate, se = std.error, p = p.value,
         ci_low = conf.low, ci_high = conf.high) |>
  arrange(p)
print(summary_tbl, n = Inf)

write_csv(summary_tbl, "output/summary_table.csv")


# ============================================================================
# 9. (Optional) LAG ANALYSIS
#
# Tourism may respond to sargassum with delay (bookings made weeks ahead).
# This re-runs the regression for the main monthly series with sargassum
# lagged by 1, 2, and 3 months.
# ============================================================================

lag_results <- list()
monthly_models <- list(
  plane       = list(d = plane_m, y = "Aantal_1"),
  nationality = list(d = nat_m,   y = "Aantal_1"),
  yachts      = list(d = yacht_m, y = "AantalJachten_1"),
  aviation    = list(d = avi_m,   y = "Passagiers_1"),
  cruise      = list(d = cru_m,   y = "Aantal_1")
)

for (nm in names(monthly_models)) {
  d <- monthly_models[[nm]]$d
  y <- monthly_models[[nm]]$y
  d <- d |> arrange(month_start)
  for (L in 1:3) {
    d_lag <- d |> mutate(positive_days_lag = lag(positive_days, L)) |>
      filter(!is.na(positive_days_lag))
    if (nrow(d_lag) < 12) next
    fml <- as.formula(sprintf("%s ~ positive_days_lag + month_fac + year_trend", y))
    m <- lm(fml, data = d_lag)
    vc <- NeweyWest(m, lag = 3, prewhite = FALSE)
    ct <- coeftest(m, vcov = vc)
    row <- ct["positive_days_lag", , drop = FALSE]
    lag_results[[paste0(nm, "_lag", L)]] <- tibble(
      dataset = nm, measure = y, lag = L, n_obs = nrow(d_lag),
      coef = row[, 1], se = row[, 2], p = row[, 4]
    )
  }
}
lag_tbl <- bind_rows(lag_results) |> arrange(dataset, lag)
print(lag_tbl, n = Inf)
write_csv(lag_tbl, "output/lag_table.csv")


# ============================================================================
# 10. (Note) Things to revisit before reporting
# ============================================================================
# - COVID 2020-2022: your sargassum series starts Dec 2020 but tourism
#   was COVID-suppressed well into 2022. Either add a COVID dummy or
#   restrict to post-recovery months.
# - Landing location: treating the whole coastline as one masks the fact
#   that a beaching at Lac Bay (tourist-dense) hits the industry harder
#   than one at the empty east shore. A spatial split of the FAI series
#   weighted by per-beach tourism volume is the proper extension.
# - Aruba / Barbados comparison: see the README for data sources.
#   A difference-in-differences with Aruba as the never-treated control
#   is the natural next analytical step.

cat("\nDone. Key outputs in ./output/:\n",
    "  sargassum_monthly.csv\n",
    "  summary_table.csv\n",
    "  lag_table.csv\n",
    "  figures/*.png\n", sep = "")
