# ============================================================================
# eda_cbs.R  --  Data mining / exploratory analysis of CBS tourism data
#
# Six tables sit alongside this script as raw CSVs (Dutch). This script
# walks through each in turn -- describing structure, plotting trends and
# seasonality, surfacing COVID effects, year-over-year growth, etc. -- and
# then combines them into a Bonaire tourism dashboard.
#
# Run section by section in RStudio (Ctrl/Cmd+Enter line by line).
# Every intermediate object stays in the global env so you can View() it.
# ============================================================================

# ---- 0. Setup ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(ggplot2)
library(scales)
library(stringr)

dir.create("eda",         showWarnings = FALSE, recursive = TRUE)
dir.create("eda/figures", showWarnings = FALSE, recursive = TRUE)
theme_set(theme_minimal(base_size = 11))

# Small helper to save a plot with a consistent filename
save_plot <- function(p, name, w = 9, h = 4) {
  ggsave(file.path("eda/figures", paste0(name, ".png")), p,
         width = w, height = h, dpi = 120)
  invisible(p)
}


# ============================================================================
# 1. Quick survey of all six tables
# ============================================================================
# What's in each file? Shapes, frequencies present, date range, NA count.

files <- c(
  aviation          = "cbs_raw/82332NED.csv",   # Caribisch NL: luchtvaart, maandcijfers luchthavens
  plane_tourism     = "cbs_raw/83104NED.csv",   # inkomend toerisme per vliegtuig
  nationality       = "cbs_raw/83191NED.csv",   # inkomend toerisme per vliegtuig; nationaliteit
  value_per_sector  = "cbs_raw/84769NED.csv",   # Bonaire; bruto toegevoegde waarde per bedrijfstak
  cruise            = "cbs_raw/85007NED.csv",   # Bonaire; cruisepassagiers
  yachts            = "cbs_raw/85015NED.csv"    # Bonaire; jachten, jachtpassagiers, ligdagen
)

survey <- lapply(files, function(f) {
  d <- read_csv(f, show_col_types = FALSE)
  freq_counts <- if ("Perioden_freq" %in% names(d))
                   paste(names(table(d$Perioden_freq)), table(d$Perioden_freq),
                         sep = "=", collapse = ", ") else "n/a"
  tibble(rows = nrow(d), cols = ncol(d), freqs = freq_counts,
         date_min = min(as.Date(d$Perioden_Date), na.rm = TRUE),
         date_max = max(as.Date(d$Perioden_Date), na.rm = TRUE))
}) |> bind_rows(.id = "table")
print(survey)
# Note: many tables contain Monthly + Quarterly + Yearly rollups stacked
# in the same file. We'll filter to the granularity we want in each section.

# ============================================================================
# 2. AVIATION  82332NED   -- monthly, by airport
# ============================================================================
# Vliegtuigbewegingen_1     aircraft movements
# TotaalPassagiersvervoer_2 total passenger transport
# AangekomenPassagiers_3    arrived passengers
# VertrokkenPassagiers_4    departed passengers
# Airports: Bonaire (Flamingo), Saba (Juancho), St. Eustatius (F.D. Roosevelt),
# plus a "Totaal Caribisch Nederland" rollup.

avi <- read_csv("cbs_raw/82332NED.csv", show_col_types = FALSE)
glimpse(avi)
unique(avi$LuchthavensCaribischNederland_label)
table(avi$Perioden_freq)  # M, Q, Y -- we want M

avi_m <- avi |>
  filter(Perioden_freq == "M") |>
  mutate(date = as.Date(Perioden_Date),
         year = year(date), month = month(date),
         airport = LuchthavensCaribischNederland_label)

# --- 2a. Long-run monthly passengers (arrivals) for each airport ----------
p_avi_arr <- avi_m |>
  filter(airport != "Totaal Caribisch Nederland") |>
  ggplot(aes(date, AangekomenPassagiers_3, colour = airport)) +
  geom_line() +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Arriving passengers by airport (monthly, 82332NED)",
       x = NULL, y = "passengers", colour = NULL)
print(p_avi_arr); save_plot(p_avi_arr, "02a_aviation_arrivals_by_airport")

# --- 2b. Bonaire arrivals -- year-month heatmap (COVID dip very visible) --
bonaire_avi <- avi_m |> filter(airport == "Flamingo Airport (Bonaire)")
p_heat <- bonaire_avi |>
  ggplot(aes(month, year, fill = AangekomenPassagiers_3)) +
  geom_tile(colour = "white") +
  scale_x_continuous(breaks = 1:12,
                     labels = c("Jan","Feb","Mar","Apr","May","Jun",
                                "Jul","Aug","Sep","Oct","Nov","Dec")) +
  scale_y_continuous(breaks = unique(bonaire_avi$year)) +
  scale_fill_viridis_c(labels = label_comma()) +
  labs(title = "Bonaire arriving passengers -- year × month heatmap",
       x = NULL, y = NULL, fill = "passengers")
print(p_heat); save_plot(p_heat, "02b_aviation_bonaire_heatmap", h = 5)

# --- 2c. Arrivals vs departures: are they balanced? -----------------------
p_balance <- bonaire_avi |>
  ggplot(aes(AangekomenPassagiers_3, VertrokkenPassagiers_4, colour = year)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.7) +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  scale_colour_viridis_c() +
  labs(title = "Bonaire arrivals vs departures (each point = one month)",
       x = "arriving", y = "departing")
print(p_balance); save_plot(p_balance, "02c_aviation_balance", w = 6, h = 5)

# --- 2d. Annual seasonal profile ------------------------------------------
season_avi <- bonaire_avi |>
  group_by(month) |>
  summarise(mean_arr = mean(AangekomenPassagiers_3, na.rm = TRUE),
            sd_arr = sd(AangekomenPassagiers_3, na.rm = TRUE),
            .groups = "drop")
p_season_avi <- ggplot(season_avi, aes(month, mean_arr)) +
  geom_ribbon(aes(ymin = mean_arr - sd_arr, ymax = mean_arr + sd_arr),
              fill = "lightblue", alpha = 0.4) +
  geom_line(linewidth = 1, colour = "#2E5C8A") +
  geom_point(colour = "#2E5C8A") +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Bonaire arriving passengers -- seasonal mean ± SD",
       x = "month", y = "mean passengers")
print(p_season_avi); save_plot(p_season_avi, "02d_aviation_seasonality", w = 7, h = 4)

# --- 2e. Annual totals & YoY growth ---------------------------------------
yearly_avi <- bonaire_avi |>
  group_by(year) |>
  summarise(total_arr = sum(AangekomenPassagiers_3, na.rm = TRUE),
            n_months  = n(), .groups = "drop") |>
  filter(n_months == 12) |>   # drop incomplete years
  mutate(yoy = (total_arr / lag(total_arr) - 1) * 100)
print(yearly_avi)
p_yoy_avi <- ggplot(yearly_avi, aes(year, yoy, fill = yoy > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#3CB371", "FALSE" = "#CD5C5C"), guide = "none") +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  scale_x_continuous(breaks = yearly_avi$year) +
  labs(title = "Bonaire arrivals -- year-over-year change", x = NULL, y = "YoY %")
print(p_yoy_avi); save_plot(p_yoy_avi, "02e_aviation_yoy", w = 7, h = 4)


# ============================================================================
# 3. TOURISM BY PLANE  83104NED  -- monthly, in thousands, by island
# ============================================================================
# Single measure: InkomendToerismePerVliegtuig_1 (x 1 000 tourists/month)
# Islands: Bonaire, Saba, Sint-Eustatius

plane <- read_csv("cbs_raw/83104NED.csv", show_col_types = FALSE)
glimpse(plane)
table(plane$CaribischNederland_label, plane$Perioden_freq)

plane_m <- plane |>
  filter(Perioden_freq == "M") |>
  mutate(date = as.Date(Perioden_Date),
         year = year(date), month = month(date),
         tourists = InkomendToerismePerVliegtuig_1 * 1000)  # was in thousands

# --- 3a. Compare three islands over time ----------------------------------
p_plane <- ggplot(plane_m, aes(date, tourists, colour = CaribischNederland_label)) +
  geom_line() +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Tourism by plane -- monthly, three islands (83104NED)",
       y = "tourists", x = NULL, colour = NULL)
print(p_plane); save_plot(p_plane, "03a_plane_three_islands")

# --- 3b. Bonaire only -- compare with aviation arrivals (sanity check) ----
bonaire_plane <- plane_m |> filter(CaribischNederland_label == "Bonaire")
join_check <- bonaire_plane |>
  select(date, tourists_plane = tourists) |>
  inner_join(bonaire_avi |> select(date, arriving_pax = AangekomenPassagiers_3), by = "date")
p_sanity <- ggplot(join_check, aes(arriving_pax, tourists_plane, colour = year(date))) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.6) +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  scale_colour_viridis_c(name = "year") +
  labs(title = "Tourism by plane vs total arrivals (both monthly, Bonaire)",
       x = "AangekomenPassagiers (82332NED)",
       y = "InkomendToerisme (83104NED)")
print(p_sanity); save_plot(p_sanity, "03b_plane_vs_aviation_sanity", w = 6, h = 5)
# Note: tourism counts a subset of passengers (residents excluded). Slope < 1
# is expected. The relationship should be roughly linear if both are sensible.

# --- 3c. Bonaire monthly trend with COVID highlight -----------------------
p_bonaire_plane <- ggplot(bonaire_plane, aes(date, tourists)) +
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2022-01-01"),
           ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.08) +
  annotate("text", x = as.Date("2021-02-01"), y = max(bonaire_plane$tourists)*0.9,
           label = "COVID", colour = "red") +
  geom_line(colour = "#2E5C8A", linewidth = 0.9) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Bonaire incoming tourism (plane), monthly",
       x = NULL, y = "tourists / month")
print(p_bonaire_plane); save_plot(p_bonaire_plane, "03c_plane_bonaire_trend")


# ============================================================================
# 4. TOURISM BY NATIONALITY  83191NED  -- YEARLY only, in PERCENT
# ============================================================================
# IMPORTANT: this table is yearly only. The unit is "% of visitors", NOT
# absolute counts. It also has a hierarchy: "Nederlands (totaal)" =
# parent of "Nederlands (Europa)" + "Nederlands (Aruba, Curaçao, St. Maarten)"
# so naive summing double-counts. Treat with care.

nat <- read_csv("cbs_raw/83191NED.csv", show_col_types = FALSE)
glimpse(nat)
table(nat$Perioden_freq)                                  # Y only
table(nat$CaribischNederland_label)                        # 3 islands
sort(unique(nat$Nationaliteit_label))                       # 15 nationality codes

nat_bo <- nat |>
  filter(CaribischNederland_label == "Bonaire") |>
  mutate(year = year(as.Date(Perioden_Date)),
         share_pct = InkomendToerismePerVliegtuig_1) |>
  select(year, nationality = Nationaliteit_label, share_pct)
glimpse(nat_bo)

# --- 4a. Sanity check: do these sum to ~100? ------------------------------
# (Spoiler: no -- 'Nederlands (totaal)' is a parent category, and individual
# nationalities are integer-rounded shares. Sum is ~150 because totaal is
# double-counted. We can drop the 'totaal' to get a cleaner ~100% picture.)
nat_bo |> group_by(year) |>
  summarise(sum_all = sum(share_pct, na.rm = TRUE),
            sum_no_totaal = sum(share_pct[nationality != "Nederlands (totaal)"], na.rm = TRUE))

# --- 4b. Top nationalities over time --------------------------------------
nat_leaf <- nat_bo |>
  filter(nationality != "Nederlands (totaal)")  # drop parent to avoid double-count
top_nats <- nat_leaf |>
  group_by(nationality) |>
  summarise(avg_share = mean(share_pct, na.rm = TRUE)) |>
  arrange(desc(avg_share))
print(top_nats)

p_nat_stack <- nat_leaf |>
  mutate(nationality = factor(nationality, levels = top_nats$nationality)) |>
  ggplot(aes(year, share_pct, fill = nationality)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d(option = "turbo") +
  scale_x_continuous(breaks = unique(nat_leaf$year)) +
  labs(title = "Bonaire visitor share by nationality (83191NED, %)",
       x = NULL, y = "share (%)", fill = NULL)
print(p_nat_stack); save_plot(p_nat_stack, "04b_nationality_stacked", w = 9, h = 5)

# --- 4c. Which nationality changed most? ----------------------------------
nat_change <- nat_leaf |>
  filter(year %in% c(min(year), max(year))) |>
  pivot_wider(names_from = year, values_from = share_pct,
              names_prefix = "y") |>
  mutate(change_pp = .data[[paste0("y", max(nat_leaf$year))]] -
                     .data[[paste0("y", min(nat_leaf$year))]]) |>
  arrange(change_pp)
print(nat_change)
# Interesting candidates: Venezolaans (collapse), Amerikaans, Dutch sub-categories.

p_nat_change <- ggplot(nat_change, aes(reorder(nationality, change_pp), change_pp,
                                        fill = change_pp > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#3CB371", "FALSE" = "#CD5C5C"), guide = "none") +
  labs(title = sprintf("Change in visitor share, %d -> %d (Bonaire)",
                        min(nat_leaf$year), max(nat_leaf$year)),
       x = NULL, y = "percentage-point change")
print(p_nat_change); save_plot(p_nat_change, "04c_nationality_change", w = 8, h = 5)


# ============================================================================
# 5. VALUE PER SECTOR  84769NED  -- YEARLY, Bonaire-only, mln USD
# ============================================================================
# Four numeric columns:
#   WaardeInWerkelijkePrijzen_1   current-price value (nominal)        mln USD
#   WaardePrijsniveau2017_2       value at 2017 prices (real)          mln USD
#   Waardemutatie_3               value change (%)
#   Volumemutatie_4               volume change (%)
# Sector "A-U Alle economische activiteiten" is the all-sector total.

vps <- read_csv("cbs_raw/84769NED.csv", show_col_types = FALSE)
glimpse(vps)
unique(vps$BedrijfstakkenBranchesSBI2008_label)

vps_clean <- vps |>
  mutate(year = year(as.Date(Perioden_Date)),
         sector = BedrijfstakkenBranchesSBI2008_label,
         nominal_mUSD = WaardeInWerkelijkePrijzen_1,
         real_mUSD    = WaardePrijsniveau2017_2,
         value_mut_pct  = Waardemutatie_3,
         volume_mut_pct = Volumemutatie_4) |>
  select(year, sector, nominal_mUSD, real_mUSD, value_mut_pct, volume_mut_pct)

# --- 5a. Total economy size over time -------------------------------------
total_econ <- vps_clean |> filter(sector == "A-U Alle economische activiteiten")
print(total_econ)
p_econ <- total_econ |>
  pivot_longer(c(nominal_mUSD, real_mUSD), names_to = "price", values_to = "value") |>
  ggplot(aes(year, value, colour = price)) +
  geom_line(linewidth = 1) + geom_point() +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Bonaire economy size (84769NED)",
       y = "mln USD", x = NULL, colour = NULL)
print(p_econ); save_plot(p_econ, "05a_econ_total", w = 7, h = 4)

# --- 5b. Sector composition stacked area ----------------------------------
sectors_only <- vps_clean |>
  filter(sector != "A-U Alle economische activiteiten") |>
  mutate(sector = str_remove(sector, "^[A-Z]+[-+]?[A-Z]?\\s+"))  # drop "A ", "D-E ", "H+J " etc.

p_sectors <- ggplot(sectors_only, aes(year, nominal_mUSD, fill = sector)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d(option = "turbo") +
  scale_x_continuous(breaks = unique(sectors_only$year)) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Bonaire value-added by sector (nominal, mln USD)",
       y = "mln USD", x = NULL, fill = NULL)
print(p_sectors); save_plot(p_sectors, "05b_sector_composition", w = 10, h = 6)

# --- 5c. Sector volume change in 2020 (COVID year) ------------------------
covid_hit <- sectors_only |> filter(year == 2020) |>
  arrange(volume_mut_pct)
print(covid_hit |> select(sector, volume_mut_pct))
p_covid <- ggplot(covid_hit, aes(reorder(sector, volume_mut_pct),
                                  volume_mut_pct, fill = volume_mut_pct > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#3CB371", "FALSE" = "#CD5C5C"), guide = "none") +
  labs(title = "Sector volume change in 2020 (% YoY) -- COVID shock",
       x = NULL, y = "volume change (%)")
print(p_covid); save_plot(p_covid, "05c_sector_covid_hit", w = 8, h = 5)

# --- 5d. Horeca (hospitality) specifically -- closest to tourism ----------
horeca <- vps_clean |> filter(sector == "I Horeca")
print(horeca)
p_horeca <- horeca |>
  pivot_longer(c(nominal_mUSD, real_mUSD), names_to = "price", values_to = "value") |>
  ggplot(aes(year, value, colour = price)) +
  geom_line(linewidth = 1) + geom_point() +
  labs(title = "Horeca (hotels & restaurants) value-added",
       y = "mln USD", x = NULL, colour = NULL)
print(p_horeca); save_plot(p_horeca, "05d_horeca", w = 7, h = 4)


# ============================================================================
# 6. CRUISE PASSENGERS  85007NED  -- monthly, Bonaire-only, in thousands
# ============================================================================

cru <- read_csv("cbs_raw/85007NED.csv", show_col_types = FALSE)
glimpse(cru)
table(cru$Perioden_freq)

cru_m <- cru |>
  filter(Perioden_freq == "M") |>
  mutate(date = as.Date(Perioden_Date),
         year = year(date), month = month(date),
         passengers = Cruisepassagiers_1 * 1000)  # was in thousands

# --- 6a. Long-run trend ---------------------------------------------------
p_cruise <- ggplot(cru_m, aes(date, passengers)) +
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2022-01-01"),
           ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.08) +
  geom_line(colour = "#7B3F9C", linewidth = 0.7) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Cruise passengers, monthly (85007NED)",
       x = NULL, y = "passengers")
print(p_cruise); save_plot(p_cruise, "06a_cruise_trend")

# --- 6b. Seasonality ------------------------------------------------------
season_cru <- cru_m |>
  filter(year < 2020 | year >= 2023) |>      # exclude COVID dip for cleaner pattern
  group_by(month) |>
  summarise(mean_pax = mean(passengers, na.rm = TRUE),
            sd_pax   = sd(passengers, na.rm = TRUE))
p_cru_season <- ggplot(season_cru, aes(month, mean_pax)) +
  geom_ribbon(aes(ymin = mean_pax - sd_pax, ymax = mean_pax + sd_pax),
              fill = "#DABFEF", alpha = 0.5) +
  geom_line(linewidth = 1, colour = "#7B3F9C") +
  geom_point(colour = "#7B3F9C") +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Cruise passengers -- seasonal mean ± SD (excl. 2020-2022)",
       x = "month", y = "mean passengers")
print(p_cru_season); save_plot(p_cru_season, "06b_cruise_seasonality", w = 7, h = 4)

# --- 6c. Annual totals ----------------------------------------------------
yearly_cru <- cru_m |>
  group_by(year) |>
  summarise(total_pax = sum(passengers, na.rm = TRUE), n_months = n()) |>
  filter(n_months == 12)
print(yearly_cru)


# ============================================================================
# 7. YACHTS  85015NED  -- monthly, Bonaire-only
# ============================================================================
# Jachtpassagiers_1          yacht passengers
# Jachten_2                  number of yachts
# GemiddeldAantalLigdagen_3  mean mooring days per yacht

yacht <- read_csv("cbs_raw/85015NED.csv", show_col_types = FALSE)
glimpse(yacht)

yacht_m <- yacht |>
  filter(Perioden_freq == "M") |>
  mutate(date = as.Date(Perioden_Date),
         year = year(date), month = month(date))

# --- 7a. Three series, same plot, scaled with secondary axis --------------
yacht_long <- yacht_m |>
  pivot_longer(c(Jachtpassagiers_1, Jachten_2, GemiddeldAantalLigdagen_3),
               names_to = "series", values_to = "value")
p_yacht <- ggplot(yacht_long, aes(date, value)) +
  geom_line(colour = "#117A65") +
  facet_wrap(~ series, scales = "free_y", ncol = 1) +
  labs(title = "Yachts: passengers, # yachts, mean mooring days (85015NED)",
       x = NULL, y = NULL)
print(p_yacht); save_plot(p_yacht, "07a_yachts_three_series", h = 7)

# --- 7b. Are mooring days inversely related to number of yachts? ----------
# (High season -> many short visits; low season -> few long-stayers?)
p_yacht_mix <- ggplot(yacht_m, aes(Jachten_2, GemiddeldAantalLigdagen_3, colour = month)) +
  geom_point(alpha = 0.7) +
  scale_colour_viridis_c(option = "C", breaks = 1:12) +
  labs(title = "Yacht count vs mean mooring days (month-coloured)",
       x = "# yachts in port", y = "mean mooring days")
print(p_yacht_mix); save_plot(p_yacht_mix, "07b_yachts_count_vs_stay", w = 6, h = 5)

# --- 7c. Seasonality of yacht passengers ----------------------------------
season_yacht <- yacht_m |>
  filter(year < 2020 | year >= 2023) |>
  group_by(month) |>
  summarise(mean_pax = mean(Jachtpassagiers_1, na.rm = TRUE))
p_yacht_season <- ggplot(season_yacht, aes(month, mean_pax)) +
  geom_col(fill = "#117A65") +
  scale_x_continuous(breaks = 1:12) +
  labs(title = "Yacht passengers -- average by month-of-year (excl. COVID)",
       x = "month", y = "mean passengers")
print(p_yacht_season); save_plot(p_yacht_season, "07c_yachts_seasonality", w = 7, h = 4)


# ============================================================================
# 8. JOINT VIEW: Bonaire tourism dashboard (monthly, 2012+)
# ============================================================================
# Build a single tibble with all monthly Bonaire tourism series, then overlay
# them on one chart (normalised so they share the y-axis).

monthly_dash <- bonaire_avi |>
  transmute(date, month_start = floor_date(date, "month"),
            aviation_arrivals = AangekomenPassagiers_3) |>
  full_join(bonaire_plane |>
              transmute(month_start = floor_date(date, "month"),
                        plane_tourism = tourists), by = "month_start") |>
  full_join(cru_m |>
              transmute(month_start = floor_date(date, "month"),
                        cruise_passengers = passengers), by = "month_start") |>
  full_join(yacht_m |>
              transmute(month_start = floor_date(date, "month"),
                        yacht_passengers = Jachtpassagiers_1), by = "month_start") |>
  select(-date) |> arrange(month_start)
glimpse(monthly_dash)

# Normalise each series to its own pre-COVID (2017-2019) mean = 100 for visual comparability
pre_covid_mean <- monthly_dash |>
  filter(month_start >= "2017-01-01" & month_start < "2020-01-01") |>
  summarise(across(-month_start, ~ mean(.x, na.rm = TRUE)))
print(pre_covid_mean)

dash_normed <- monthly_dash |>
  mutate(across(-month_start,
                ~ .x / pre_covid_mean[[cur_column()]] * 100))
dash_long <- dash_normed |>
  pivot_longer(-month_start, names_to = "series", values_to = "index")

p_dash <- ggplot(dash_long, aes(month_start, index, colour = series)) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey50") +
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2022-01-01"),
           ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.06) +
  geom_line(alpha = 0.85) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Bonaire monthly tourism series -- indexed (2017-2019 mean = 100)",
       x = NULL, y = "index", colour = NULL)
print(p_dash); save_plot(p_dash, "08_bonaire_dashboard_indexed", w = 11, h = 5)


# ============================================================================
# 9. Recovery from COVID: who came back fastest?
# ============================================================================
# Compare the most recent complete year to 2019 across all four series.

# A year is "complete" if all 12 months are present.
complete_year <- function(df, date_col = "date") {
  df |> mutate(yr = year(.data[[date_col]])) |>
    count(yr) |> filter(n == 12) |> pull(yr)
}
yrs_avi   <- complete_year(bonaire_avi)
yrs_cru   <- complete_year(cru_m)
yrs_yacht <- complete_year(yacht_m)
yrs_plane <- complete_year(bonaire_plane)
last_full <- min(max(yrs_avi),   max(yrs_cru),
                 max(yrs_yacht), max(yrs_plane))
cat("Last full year (all series have 12 months):", last_full, "\n")

recovery <- bind_rows(
  bonaire_avi |> group_by(year) |>
    summarise(total = sum(AangekomenPassagiers_3, na.rm = TRUE), .groups = "drop") |>
    mutate(series = "aviation_arrivals"),
  cru_m |> group_by(year) |>
    summarise(total = sum(passengers, na.rm = TRUE), .groups = "drop") |>
    mutate(series = "cruise"),
  yacht_m |> group_by(year) |>
    summarise(total = sum(Jachtpassagiers_1, na.rm = TRUE), .groups = "drop") |>
    mutate(series = "yachts"),
  bonaire_plane |> group_by(year) |>
    summarise(total = sum(tourists, na.rm = TRUE), .groups = "drop") |>
    mutate(series = "plane_tourism")
) |>
  filter(year %in% c(2019, last_full)) |>
  pivot_wider(names_from = year, values_from = total, names_prefix = "y") |>
  mutate(recovery_pct = .data[[paste0("y", last_full)]] /
                       .data[["y2019"]] * 100)
print(recovery)

p_recovery <- ggplot(recovery, aes(reorder(series, recovery_pct), recovery_pct)) +
  geom_col(fill = "#2E5C8A") +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "red") +
  geom_text(aes(label = sprintf("%.0f%%", recovery_pct)),
            hjust = -0.1, size = 4) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(recovery$recovery_pct) * 1.1)) +
  labs(title = sprintf("Bonaire tourism recovery: %d vs 2019", last_full),
       x = NULL, y = "% of 2019 level")
print(p_recovery); save_plot(p_recovery, "09_recovery_2019", w = 8, h = 4)


# ============================================================================
# 10. Year × month heatmaps for all monthly Bonaire series, on one canvas
# ============================================================================
# Spotting structural breaks visually.

dash_heat <- monthly_dash |>
  pivot_longer(-month_start, names_to = "series", values_to = "value") |>
  mutate(year = year(month_start), month = month(month_start)) |>
  group_by(series) |>
  # Normalize 0-1 within each series so the four panels are visually comparable
  mutate(value_norm = (value - min(value, na.rm = TRUE)) /
                       (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))) |>
  ungroup()

p_heat_all <- ggplot(dash_heat, aes(month, year, fill = value_norm)) +
  geom_tile(colour = "white", linewidth = 0.1) +
  facet_wrap(~ series, ncol = 2) +
  scale_x_continuous(breaks = 1:12, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(2012, 2026, 2), trans = "reverse") +
  scale_fill_viridis_c(na.value = "grey90", labels = label_percent(),
                       name = "% of series\nmax") +
  labs(title = "Year-month heatmaps -- Bonaire monthly tourism (each series normalised)",
       x = NULL, y = NULL)
print(p_heat_all); save_plot(p_heat_all, "10_heatmaps_all", w = 11, h = 7)


# ============================================================================
# 11. Save key cleaned tables for downstream use
# ============================================================================
write_csv(monthly_dash,  "eda/bonaire_monthly_dashboard.csv")
write_csv(yearly_avi,    "eda/aviation_yearly.csv")
write_csv(yearly_cru,    "eda/cruise_yearly.csv")
write_csv(nat_leaf,      "eda/nationality_yearly_bonaire.csv")
write_csv(vps_clean,     "eda/value_per_sector_bonaire.csv")
write_csv(recovery,      "eda/recovery_table.csv")

cat("\nDone. Open eda/figures/ for plots and eda/*.csv for cleaned tables.\n")
