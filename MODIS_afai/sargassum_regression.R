# =============================================================================
# Sargassum Influx & Tourism Regression — Bonaire (+ Aruba control)
# Data sources:
#   - Sargassum proxy : NOAA/AOML ERDDAP — USF AFAI 1-day fields (2016–present)
#   - Bonaire tourism : Statistics Netherlands (CBS) via cbsodataR
#   - Aruba tourism   : World Bank API via WDI (annual) /
#                       CBS Aruba manual download (monthly, see note below)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. PACKAGES
# -----------------------------------------------------------------------------
pkgs <- c("cbsodataR", "WDI", "httr", "dplyr", "tidyr", "lubridate",
          "ggplot2", "patchwork", "lmtest", "sandwich", "scales", "readr")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------
# Bounding boxes (lon_min, lat_min, lon_max, lat_max)
# Bonaire: focus on eastern coast / Lac Bay area per van der Geest (2024)
BONAIRE_LAT <- c(12.00, 12.35)
BONAIRE_LON <- c(-68.40, -68.10)

# Aruba: west/south coasts face the Caribbean (primary tourist beaches)
ARUBA_LAT   <- c(12.35, 12.65)
ARUBA_LON   <- c(-70.10, -69.85)

# AFAI positive threshold — values above this indicate floating algae
# Wang & Hu (2016) suggest 0.0 as minimum; 0.0005 reduces sun-glint noise
AFAI_THRESHOLD <- 0.0005

# Date range — AFAI dataset starts 2016-06-18
AFAI_START <- "2016-07-01"   # first full month available
AFAI_END   <- format(Sys.Date() - 1, "%Y-%m-%d")

# -----------------------------------------------------------------------------
# 2. ERDDAP AFAI DOWNLOAD
# Endpoint: NOAA CoastWatch Caribbean / AOML
# Dataset : noaa_aoml_atlantic_oceanwatch_AFAI_1D
# Docs    : https://cwcgom.aoml.noaa.gov/erddap/griddap/
#           noaa_aoml_atlantic_oceanwatch_AFAI_1D.html
# -----------------------------------------------------------------------------
ERDDAP_BASE <- paste0(
  "https://cwcgom.aoml.noaa.gov/erddap/griddap/",
  "noaa_aoml_atlantic_oceanwatch_AFAI_1D.csv"
)

#' Download daily AFAI for a bounding box, return tidy data frame.
#' Downloads year-by-year to keep individual requests small (~3 MB each).
#'
#' @param lat_range  numeric(2) c(lat_min, lat_max)
#' @param lon_range  numeric(2) c(lon_min, lon_max)
#' @param start_date character "YYYY-MM-DD"
#' @param end_date   character "YYYY-MM-DD"
#' @return data.frame with columns: date, afai_mean, afai_positive_frac
download_afai <- function(lat_range, lon_range, start_date, end_date) {
  years <- seq(year(start_date), year(end_date))
  all_rows <- vector("list", length(years))

  for (i in seq_along(years)) {
    yr    <- years[i]
    t_min <- if (yr == year(start_date)) start_date else paste0(yr, "-01-01")
    t_max <- if (yr == year(end_date))   end_date   else paste0(yr, "-12-31")

    # ERDDAP griddap URL — coordinate constraints use parentheses notation:
    # variable[(t_start):(t_end)][(lat_min):(lat_max)][(lon_min):(lon_max)]
    query <- sprintf(
      "AFAI[(%sT12:00:00Z):(%sT12:00:00Z)][(%s):(%s)][(%s):(%s)]",
      t_min, t_max,
      lat_range[1], lat_range[2],
      lon_range[1],  lon_range[2]
    )
    url <- paste0(ERDDAP_BASE, "?", query)

    message(sprintf("  Downloading AFAI for bbox %.2f–%.2f°N, %.2f–%.2f°E, year %d ...",
                    lat_range[1], lat_range[2], lon_range[1], lon_range[2], yr))

    resp <- tryCatch(
      GET(url, timeout(120)),
      error = function(e) { message("  Request failed: ", e$message); NULL }
    )

    if (is.null(resp) || status_code(resp) != 200) {
      message("  HTTP ", status_code(resp), " — skipping year ", yr)
      next
    }

    # ERDDAP returns 2 header rows (names + units) — skip row 2
    raw_text <- content(resp, as = "text", encoding = "UTF-8")
    df <- tryCatch({
      lines <- strsplit(raw_text, "\n")[[1]]
      lines <- lines[c(1, 3:length(lines))]          # drop units row
      read_csv(paste(lines, collapse = "\n"),
               col_types = cols(.default = col_double(),
                                time = col_character()),
               show_col_types = FALSE)
    }, error = function(e) {
      message("  Parse error: ", e$message); NULL
    })

    if (is.null(df) || nrow(df) == 0) next

    colnames(df) <- tolower(colnames(df))             # normalise column names

    all_rows[[i]] <- df |>
      mutate(date = as.Date(substr(time, 1, 10))) |>
      select(date, afai = afai)
  }

  bind_rows(all_rows)
}

#' Aggregate daily pixel-level AFAI to a monthly sargassum index per island.
#'
#' Monthly metrics returned:
#'   afai_mean         — mean AFAI over all valid (non-NA) pixels & days
#'   afai_pos_mean     — mean AFAI for pixels > threshold (sargassum pixels only)
#'   sarg_pixel_frac   — fraction of pixel-days with AFAI > threshold
#'   n_days            — days with data in that month
monthly_afai <- function(daily_df, threshold = AFAI_THRESHOLD) {
  daily_df |>
    filter(!is.na(afai)) |>
    mutate(month = floor_date(date, "month"),
           is_sarg = afai > threshold) |>
    group_by(month) |>
    summarise(
      afai_mean       = mean(afai,             na.rm = TRUE),
      afai_pos_mean   = mean(afai[is_sarg],    na.rm = TRUE),
      sarg_pixel_frac = mean(is_sarg,          na.rm = TRUE),
      n_days          = n_distinct(date),
      .groups = "drop"
    ) |>
    # Require at least 4 days with data; AFAI saturates under clouds
    filter(n_days >= 4)
}

# -----------------------------------------------------------------------------
# 3. DOWNLOAD & PROCESS AFAI
# -----------------------------------------------------------------------------
message("=== Downloading AFAI for Bonaire ===")
bonaire_daily <- download_afai(BONAIRE_LAT, BONAIRE_LON, AFAI_START, AFAI_END)
bonaire_afai  <- monthly_afai(bonaire_daily)

message("=== Downloading AFAI for Aruba ===")
aruba_daily   <- download_afai(ARUBA_LAT, ARUBA_LON, AFAI_START, AFAI_END)
aruba_afai    <- monthly_afai(aruba_daily)

message(sprintf("Bonaire: %d monthly AFAI records", nrow(bonaire_afai)))
message(sprintf("Aruba  : %d monthly AFAI records", nrow(aruba_afai)))

# -----------------------------------------------------------------------------
# 4. CBS — BONAIRE TOURISM DATA
# Using cbsodataR's cbs_add_date_column() to handle CBS period codes cleanly.
# CBS period code format: "2023MM04" = April 2023, "2023JJ00" = full year 2023
# -----------------------------------------------------------------------------
message("=== Loading CBS tourism data for Bonaire ===")

#' Load and clean a CBS Caribbean NL tourism table.
#' Returns monthly data for one island.
load_cbs_tourism <- function(table_id, island_pattern = "Bonaire") {
  meta <- cbs_get_meta(table_id)

  # Print DataProperties so you can identify the right value column
  message("  Table columns:")
  print(meta$DataProperties[, c("Key", "Title", "Unit")])

  raw <- cbs_get_data(table_id)

  # cbs_add_date_column converts "2023MM04" -> Date, adds 'Date' and 'freq' columns
  df <- cbs_add_date_column(raw, date_type = "Date")

  # Identify island/region column (CBS uses Dutch: 'Eilanden', English: 'Islands' etc.)
  region_col <- names(df)[sapply(names(df), function(n)
    any(c("eiland", "island", "territory", "caribbean", "gebied") %in% tolower(n)))]

  if (length(region_col) > 0) {
    df <- df |> filter(grepl(island_pattern, .data[[region_col[1]]], ignore.case = TRUE))
    message(sprintf("  Filtered to '%s' using column '%s': %d rows",
                    island_pattern, region_col[1], nrow(df)))
  } else {
    message("  Warning: no region column found — keeping all rows")
  }

  # Keep only monthly rows (freq == "M" after cbs_add_date_column)
  df <- df |> filter(freq == "M")

  # Identify the primary numeric value column (first non-ID numeric column)
  skip_cols <- c("ID", names(df)[sapply(names(df), function(n) inherits(df[[n]], "Date"))])
  num_cols  <- names(df)[sapply(names(df), is.numeric)]
  num_cols  <- setdiff(num_cols, c("ID"))

  if (length(num_cols) == 0) stop("No numeric columns found in table ", table_id)

  # Show candidate columns so researcher can verify the right one
  message("  Numeric candidates: ", paste(num_cols, collapse = ", "))

  # Use first numeric column; researcher can override by changing [1]
  val_col <- num_cols[1]
  message("  Using value column: ", val_col)

  df |>
    select(month = Date, visitors = all_of(val_col)) |>
    filter(!is.na(visitors), visitors > 0)
}
# A043194 = bonaire airport
bonaire_tourism <- tryCatch(
  load_cbs_tourism("82332NED", "Bonaire"),
  error = function(e) {
    message("CBS load failed: ", e$message)

    # --- FALLBACK: Search CBS for the right table ---
    message("Searching CBS for Caribbean NL tourism tables...")
    ds <- cbs_get_datasets() |>
      filter(grepl("Caribisch|bonaire|passagiers|toerisme", Title, ignore.case = TRUE))
    message("Candidate tables:")
    print(ds[, c("Identifier", "Title", "Updated")])
    return(data.frame())
  }
)

raw <- cbs_add_date_column(raw, date_type = "Date")

# A043194 = bonaire airport
bonaire_tourism <- raw %>%
  filter(LuchthavensCaribischNederland == "A043194")

# -----------------------------------------------------------------------------
# 4b. ARUBA TOURISM DATA
#
# Aruba is a constituent country of the Kingdom of Netherlands — NOT a Dutch
# special municipality — so it is NOT in the CBS Netherlands open data portal.
#
# Option A (used here): World Bank Development Indicators — annual arrivals.
#   Indicator: ST.INT.ARVL  (International inbound tourists)
#   Useful for year-level comparison with annual aggregated AFAI.
#
# Option B (monthly, manual): CBS Aruba publishes monthly tourist arrivals at
#   https://cbs.aw/wp/index.php/category/tourism/
#   Download their Excel/CSV files, read with readxl::read_excel(), and merge.
#   Column: "Number of stayover" by month.
#
# Option C (monthly, semi-automated): Aruba Tourism Authority press releases
#   https://www.aruba.com/us/media/press-releases
# -----------------------------------------------------------------------------
message("=== Loading Aruba tourism data (World Bank, annual) ===")
aruba_wdi <- tryCatch({
  WDI(country = "AW",             # ISO2 code for Aruba
      indicator = "ST.INT.ARVL",  # International tourist arrivals
      start = 2010,
      end   = year(Sys.Date())) |>
    as_tibble() |>
    rename(year = year, arrivals_aruba = ST.INT.ARVL) |>
    filter(!is.na(arrivals_aruba)) |>
    select(year, arrivals_aruba)
}, error = function(e) {
  message("WDI download failed: ", e$message)
  data.frame()
})

message(sprintf("Aruba WDI: %d annual records", nrow(aruba_wdi)))

# Aggregate monthly Aruba AFAI to annual for comparison with annual WDI data
aruba_afai_annual <- aruba_afai |>
  mutate(year = year(month)) |>
  group_by(year) |>
  summarise(afai_mean_annual = mean(afai_mean, na.rm = TRUE),
            sarg_frac_annual = mean(sarg_pixel_frac, na.rm = TRUE),
            n_months = n(), .groups = "drop") |>
  filter(n_months >= 9)   # only include years with near-complete data

# -----------------------------------------------------------------------------
# 5. MERGE DATASETS
# -----------------------------------------------------------------------------

# --- 5a. Monthly: Bonaire AFAI + Bonaire CBS tourism ---
bonaire_monthly <- if (nrow(bonaire_tourism) > 0 && nrow(bonaire_afai) > 0) {
  inner_join(
    bonaire_afai,
    bonaire_tourism |> rename(month = Perioden_Date),
    by = "month"
  ) |> arrange(month)
} else {
  message("WARNING: Cannot create Bonaire monthly merged dataset.")
  data.frame()
}

message(sprintf("Bonaire monthly merged: %d observations (%s to %s)",
                nrow(bonaire_monthly),
                if (nrow(bonaire_monthly) > 0) as.character(min(bonaire_monthly$month)) else "NA",
                if (nrow(bonaire_monthly) > 0) as.character(max(bonaire_monthly$month)) else "NA"))

# --- 5b. Annual: Aruba AFAI + Aruba WDI tourism ---
aruba_annual <- if (nrow(aruba_wdi) > 0 && nrow(aruba_afai_annual) > 0) {
  inner_join(aruba_wdi, aruba_afai_annual, by = "year")
} else {
  data.frame()
}

# -----------------------------------------------------------------------------
# 6. VISUALISATION
# -----------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40"),
        axis.title = element_text(size = 11)))

# --- 6a. AFAI time series ---
p_afai <- bind_rows(
  bonaire_afai |> mutate(island = "Bonaire"),
  aruba_afai   |> mutate(island = "Aruba")
) |>
  ggplot(aes(month, sarg_pixel_frac, color = island)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, linewidth = 1.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_color_manual(values = c(Bonaire = "#0077B6", Aruba = "#E63946")) +
  labs(title = "Monthly Sargassum Pixel Fraction — AFAI > threshold",
       subtitle = paste0("Source: NOAA/AOML USF AFAI 1D (threshold = ", AFAI_THRESHOLD, ")"),
       x = NULL, y = "Fraction of pixels with AFAI > threshold",
       color = NULL) +
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2021-12-01"),
           ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "grey30") +
  annotate("text", x = as.Date("2021-01-01"), y = Inf,
           label = "COVID-19", vjust = 1.5, size = 3, color = "grey40")

# --- 6b. Bonaire tourism time series ---
p_tourism <- if (nrow(bonaire_tourism) > 0) {
  bonaire_tourism |>
    ggplot(aes(Perioden_Date, AangekomenPassagiers_3)) +
    geom_col(fill = "#0077B6", alpha = 0.75) +
    geom_smooth(method = "loess", span = 0.2, se = FALSE,
                color = "navy", linewidth = 1) +
    scale_y_continuous(labels = comma_format()) +
    labs(title = "Bonaire Monthly Tourist Arrivals by Air",
         subtitle = "Source: Statistics Netherlands (CBS), table 83104ENG",
         x = NULL, y = "Arrivals")
} else {
  ggplot() + labs(title = "Bonaire tourism data not loaded")
}

# --- 6c. Scatter: AFAI vs visitors (Bonaire monthly) ---
p_scatter_bonaire <- if (nrow(bonaire_monthly) >= 5) {
  bonaire_monthly |>
    ggplot(aes(sarg_pixel_frac, log10(AangekomenPassagiers_3))) +
    geom_point(alpha = 0.6, color = "#0077B6", size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "navy") +
    scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
    scale_y_continuous(labels = comma_format()) +
    labs(title = "Bonaire: AFAI Fraction vs. Monthly Arrivals",
         x = "Sargassum pixel fraction", y = "Monthly arrivals")
} else {
  ggplot() + labs(title = "Insufficient data for Bonaire scatter")
}

# --- 6d. Scatter: annual Aruba (control) ---
p_scatter_aruba <- if (nrow(aruba_annual) >= 4) {
  aruba_annual |>
    ggplot(aes(sarg_frac_annual, arrivals_aruba / 1e3)) +
    geom_point(color = "#E63946", size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "darkred") +
    geom_text(aes(label = year), vjust = -0.8, size = 3, color = "grey40") +
    scale_y_continuous(labels = comma_format(suffix = "k")) +
    labs(title = "Aruba (control): Annual AFAI vs. Tourist Arrivals",
         subtitle = "Source: World Bank ST.INT.ARVL",
         x = "Mean annual sargassum fraction", y = "Annual arrivals (thousands)")
} else {
  ggplot() + labs(title = "Insufficient data for Aruba scatter")
}

# Combine and save
combined_plot <- (p_afai / p_tourism) | (p_scatter_bonaire / p_scatter_aruba)
ggsave("sargassum_tourism_overview.png", combined_plot,
       width = 16, height = 10, dpi = 150)
message("Plot saved: sargassum_tourism_overview.png")

# -----------------------------------------------------------------------------
# 7. REGRESSION ANALYSIS
# -----------------------------------------------------------------------------

#' Run OLS with HAC (Newey-West) standard errors and a lag.
#' Newey-West is appropriate here: the residuals are likely both
#' heteroskedastic and autocorrelated (monthly time series).
#'
#' @param df        data.frame with columns `sarg_pixel_frac` and `visitors`
#' @param lag_m     integer months to lag the AFAI predictor (0 = no lag)
#' @param log_y     logical whether to log-transform visitors
#' @return list(model, coeftest_hac, df_used)
run_ols <- function(df, lag_m = 0, log_y = TRUE) {
  d <- df |>
    arrange(month) |>
    mutate(
      sarg_lag   = lag(sarg_pixel_frac, lag_m),
      y          = if (log_y) log(visitors + 1) else visitors,
      # Season controls: month dummies (important — tourism is highly seasonal)
      month_num  = month(month)
    ) |>
    filter(!is.na(sarg_lag), !is.na(y))

  if (nrow(d) < 10) {
    message("  Insufficient observations (n=", nrow(d), ") for lag=", lag_m)
    return(NULL)
  }

  # OLS with month-of-year fixed effects to remove seasonality
  formula_str <- "y ~ sarg_lag + factor(month_num)"
  mod <- lm(as.formula(formula_str), data = d)

  # HAC standard errors (Newey-West, max lag = 4 months)
  hac <- coeftest(mod, vcov = NeweyWest(mod, lag = 4, prewhite = FALSE))

  list(model = mod, hac = hac, data = d, lag = lag_m, log_y = log_y)
}

# --- 7a. Bonaire monthly regression ---
if (nrow(bonaire_monthly) >= 10) {
  message("\n=== Bonaire OLS Regression (monthly, log visitors ~ AFAI + month FE) ===")

  results_bonaire <- lapply(0:3, function(lag) run_ols(bonaire_monthly, lag_m = lag))
  results_bonaire <- results_bonaire[!sapply(results_bonaire, is.null)]

  for (res in results_bonaire) {
    cat(sprintf("\n--- Lag = %d month(s) ---\n", res$lag))
    print(res$hac)
    r2 <- summary(res$model)$r.squared
    n  <- nrow(res$data)
    cat(sprintf("R² = %.3f  |  n = %d\n", r2, n))
  }

  # Extract sarg_lag coefficient and p-value across lags for a summary table
  lag_summary <- lapply(results_bonaire, function(res) {
    h <- res$hac
    tibble(
      lag      = res$lag,
      coef     = h["sarg_lag", "Estimate"],
      se_hac   = h["sarg_lag", "Std. Error"],
      p_hac    = h["sarg_lag", "Pr(>|t|)"],
      r2       = summary(res$model)$r.squared,
      n        = nrow(res$data)
    )
  }) |> bind_rows()

  cat("\n=== Bonaire lag sensitivity summary ===\n")
  print(lag_summary, n = Inf)

  # Publication-quality lag sensitivity plot
  p_lags <- lag_summary |>
    mutate(
      ci_lo = coef - 1.96 * se_hac,
      ci_hi = coef + 1.96 * se_hac,
      sig   = p_hac < 0.05
    ) |>
    ggplot(aes(factor(lag), coef, color = sig)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), size = 0.8, linewidth = 1) +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#0077B6"),
                       labels = c("n.s.", "p < 0.05"),
                       name = NULL) +
    labs(title = "Bonaire: AFAI coefficient on log(arrivals) by lag",
         subtitle = "HAC (Newey-West) 95% CI | Controls: month-of-year FE",
         x = "Lag (months)", y = "Coefficient (log scale)")

  ggsave("bonaire_lag_sensitivity.png", p_lags,
         width = 7, height = 5, dpi = 150)
  message("Lag sensitivity plot saved: bonaire_lag_sensitivity.png")

} else {
  message("Not enough Bonaire monthly observations for regression.")
}

# --- 7b. Aruba annual regression (control check) ---
if (nrow(aruba_annual) >= 5) {
  message("\n=== Aruba annual regression (control: expect weak/no sargassum effect) ===")
  aruba_mod <- lm(log(arrivals_aruba + 1) ~ sarg_frac_annual, data = aruba_annual)
  print(summary(aruba_mod))
  cat("Interpretation: If Aruba shows NO significant sargassum coefficient while\n")
  cat("Bonaire DOES, that supports a causal interpretation for Bonaire.\n")
}

# --- 7c. Difference-in-Differences framing (bonus if monthly Aruba data available) ---
# If you obtain monthly Aruba CBS/ATA data, uncomment and adapt:
#
# aruba_monthly <- read_csv("aruba_monthly_arrivals.csv") |>
#   rename(month = date, visitors_aruba = arrivals)
#
# did_df <- bonaire_monthly |>
#   left_join(aruba_monthly, by = "month") |>
#   mutate(
#     log_ratio = log((visitors + 1) / (visitors_aruba + 1)),
#     month_num = month(month)
#   )
#
# did_mod <- lm(log_ratio ~ sarg_pixel_frac + factor(month_num), data = did_df)
# print(coeftest(did_mod, vcov = NeweyWest(did_mod, lag = 4, prewhite = FALSE)))
# # Coefficient on sarg_pixel_frac = DiD estimate of sargassum on relative tourism

# -----------------------------------------------------------------------------
# 8. SAVE CLEAN DATA FOR STAKEHOLDER SLIDES
# -----------------------------------------------------------------------------
if (nrow(bonaire_monthly) > 0) {
  write_csv(bonaire_monthly, "bonaire_sargassum_tourism_monthly.csv")
  message("Saved: bonaire_sargassum_tourism_monthly.csv")
}
if (nrow(aruba_afai) > 0) {
  write_csv(aruba_afai, "aruba_afai_monthly.csv")
  message("Saved: aruba_afai_monthly.csv")
}
if (nrow(aruba_annual) > 0) {
  write_csv(aruba_annual, "aruba_annual_merged.csv")
  message("Saved: aruba_annual_merged.csv")
}

message("\n=== Done. Check the PNG plots and CSV files. ===")

# -----------------------------------------------------------------------------
# NOTES FOR RESEARCHER
# -----------------------------------------------------------------------------
# 1. AFAI LIMITATIONS:
#    - Daily AFAI saturates under clouds/sun-glint; months with < 4 clear days
#      are excluded. Cloud cover is high in Caribbean summer.
#    - AFAI detects ALL floating algae, not only sargassum. In Bonaire's small
#      bbox the dominant floating material is sargassum, but cite this caveat.
#    - The dataset starts June 2016, limiting the pre-invasion baseline.
#      CBS tourism data from 2012 can still be used as pre-period context
#      (just without AFAI — show it descriptively in slide 1).
#
# 2. ARUBA MONTHLY DATA (for full DiD):
#    - CBS Aruba: https://cbs.aw/wp/index.php/category/tourism/
#      Download "Table 1.1 Number of stayover" Excel files, combine years,
#      read into R with: readxl::read_excel("filename.xlsx", skip = 2)
#    - Aruba Tourism Authority: https://www.aruba.com/us/media/press-releases
#      Monthly press releases with stayover + cruise numbers.
#    - Once loaded, use the DiD block above (Section 7c).
#
# 3. SEASONALITY:
#    - Tourism and sargassum BOTH peak in summer — month FE controls are
#      essential to avoid spurious correlation. They are included above.
#
# 4. CBS TABLE DISCOVERY:
#    Run this to search for additional tourism tables:
#      cbs_get_datasets() |>
#        filter(grepl("tourist|arrival|caribbean|bonaire", Title, ignore.case=TRUE),
#               grepl("ENG|eng", Identifier)) |>
#        select(Identifier, Title, Updated)
