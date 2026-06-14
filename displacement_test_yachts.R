# =============================================================================
# Displacement / refuge hypothesis: do regional sargassum influxes push yachts
# towards Bonaire?
#
# Hypothesis: the positive own-sargassum coefficient on Bonaire yacht counts is
# spurious. When sargassum floods the wider eastern Caribbean, yachts relocate
# from high-exposure islands (here proxied by Barbados) towards Bonaire, so
# Bonaire yacht activity rises with REGIONAL pressure and, in particular, with
# Bonaire being relatively cleaner than its neighbours.
#
# Sargassum variables are taken ENTIRELY from the detected-event dataset
# (detected_events.csv) - the same connected-component + temporal-persistence
# events used elsewhere in the thesis. We expand each event across the calendar
# months it spans and build two monthly exposure measures per island:
#   - event_days : days within the month falling inside a detected event
#   - peak_cov   : peak coverage fraction of events overlapping the month
# (mean_cov is also carried for reference).
#
# Design: single monthly time series for Bonaire (yacht data are Bonaire-only),
# with Bonaire's own event exposure and Barbados' event exposure entered side by
# side, plus a relative-exposure gap. Month-of-year fixed effects throughout, a
# year/trend robustness layer, negative-binomial counts, Newey-West HAC errors.
# =============================================================================

library(dplyr)
library(lubridate)
library(MASS)         # glm.nb
library(sandwich)     # NeweyWest
library(lmtest)       # coeftest
library(modelsummary)

# -----------------------------------------------------------------------------
# 1. Monthly sargassum exposure from the detected-event dataset
# -----------------------------------------------------------------------------
events <- read.csv("detected_events.csv", stringsAsFactors = FALSE) %>%
  mutate(island = trimws(island),
         start  = as.Date(start),
         end    = as.Date(end))

# Expand a single event into its monthly day-contributions. Events within an
# island are sequential and non-overlapping, so monthly day counts simply sum.
expand_event <- function(isl, s, e, pk, mc) {
  m_starts <- seq(floor_date(s, "month"), floor_date(e, "month"), by = "month")
  m_ends   <- ceiling_date(m_starts, "month") - 1          # last day of each month
  data.frame(
    island   = isl,
    year     = year(m_starts),
    month    = month(m_starts),
    ev_days  = as.integer(pmin(e, m_ends) - pmax(s, m_starts)) + 1L,
    peak_cov = pk,
    mean_cov = mc,
    stringsAsFactors = FALSE
  )
}

monthly_events <- lapply(seq_len(nrow(events)), function(i) {
  expand_event(events$island[i], events$start[i], events$end[i],
               events$peak_coverage[i], events$mean_coverage[i])
}) %>% bind_rows()

# Collapse to one row per island-month
sarg <- monthly_events %>%
  group_by(island, year, month) %>%
  summarise(
    event_days = sum(ev_days),                       # days-in-event (frequency)
    peak_cov   = max(peak_cov),                       # peak intensity
    mean_cov   = sum(mean_cov * ev_days) / sum(ev_days),  # event-day-weighted mean
    n_events   = n(),
    .groups = "drop"
  )

bon_sarg <- sarg %>% filter(island == "Bonaire")  %>% dplyr::select(-island, -n_events)
bar_sarg <- sarg %>% filter(island == "Barbados") %>% dplyr::select(-island, -n_events)
names(bon_sarg)[-(1:2)] <- paste0("bon_", names(bon_sarg)[-(1:2)])
names(bar_sarg)[-(1:2)] <- paste0("bar_", names(bar_sarg)[-(1:2)])

# -----------------------------------------------------------------------------
# 2. Monthly Bonaire yacht outcome (CBS 85015NED)
# -----------------------------------------------------------------------------
# Keep only the monthly frequency rows (drop quarterly "KW" and annual totals).
yacht <- read.csv("cbs_raw/85015NED.csv", stringsAsFactors = FALSE) %>%
  filter(Perioden_freq == "M") %>%
  mutate(date  = as.Date(Perioden_Date),
         year  = year(date),
         month = month(date),
         yacht_n   = Jachten_2,
         yacht_pax = Jachtpassagiers_1,
         mooring   = GemiddeldAantalLigdagen_3) %>%
  dplyr::select(year, month, yacht_n, yacht_pax, mooring)

# -----------------------------------------------------------------------------
# 3. Build the merged monthly panel and the displacement variables
# -----------------------------------------------------------------------------
# Window: 2021-01 .. 2025-12 (60 months). Months with no detected event enter as
# genuine zeros (they were scanned by the detector and produced no event), which
# preserves the contrast between non-event and event months.

panel <- yacht %>%
  filter(year >= 2021, year <= 2025) %>%
  left_join(bon_sarg, by = c("year", "month")) %>%
  left_join(bar_sarg, by = c("year", "month")) %>%
  mutate(across(c(bon_event_days, bon_peak_cov, bon_mean_cov,
                  bar_event_days, bar_peak_cov, bar_mean_cov),
                ~ tidyr::replace_na(., 0))) %>%
  arrange(year, month) %>%
  mutate(
    # Relative-exposure gap: positive = Barbados worse than Bonaire = Bonaire is
    # the relatively cleaner refuge. The refuge story predicts a positive sign.
    gap_event_days = bar_event_days - bon_event_days,
    gap_peak_cov   = bar_peak_cov   - bon_peak_cov,
    gap_mean_cov   = bar_mean_cov   - bon_mean_cov,
    log_yacht   = log(yacht_n),
    log_mooring = log(mooring),
    month_f   = factor(month),                 # month-of-year fixed effects
    year_f    = factor(year),                  # year fixed effects (robustness)
    t         = (year - 2021) * 12 + month     # linear trend (robustness)
  )

cat("Panel:", nrow(panel), "months,", min(panel$year), "-", max(panel$year), "\n")
cat("corr(bon_event_days, bar_event_days) =",
    round(cor(panel$bon_event_days, panel$bar_event_days), 3),
    " corr(bon_peak_cov, bar_peak_cov) =",
    round(cor(panel$bon_peak_cov, panel$bar_peak_cov), 3), "\n")
cat("corr(bon_event_days, trend) =",
    round(cor(panel$bon_event_days, panel$t), 3),
    " corr(yacht_n, trend) =",
    round(cor(panel$yacht_n, panel$t), 3), "\n\n")

# -----------------------------------------------------------------------------
# 4. Helpers: main-model families reported with Newey-West HAC (lag 3)
# -----------------------------------------------------------------------------
# Yacht count is an overdispersed count -> negative binomial with a log link.
# Mooring days is a continuous outcome -> log-linear OLS (log on the LHS).
# Both report Newey-West HAC standard errors, exactly as in the main models.
nb_hac <- function(formula, data, lag = 3) {
  m  <- glm.nb(formula, data = data)
  ct <- coeftest(m, vcov. = NeweyWest(m, lag = lag, prewhite = FALSE))
  list(model = m, ct = ct)
}

lm_hac <- function(formula, data, lag = 3) {
  m  <- lm(formula, data = data)
  ct <- coeftest(m, vcov. = NeweyWest(m, lag = lag, prewhite = FALSE))
  list(model = m, ct = ct)
}

# -----------------------------------------------------------------------------
# 5. Model A - relative-exposure gap (Barbados minus Bonaire), month FE
# -----------------------------------------------------------------------------
# Refuge prediction: positive sign (Bonaire cleaner relative to Barbados ->
# more yachts seek it out). Each gap measure is run against both yacht outcomes.
# Yacht count: negative binomial (log link).
mA_cnt_days  <- nb_hac(yacht_n     ~ gap_event_days + month_f, panel)
mA_cnt_peak  <- nb_hac(yacht_n     ~ gap_peak_cov   + month_f, panel)
mA_cnt_mean  <- nb_hac(yacht_n     ~ gap_mean_cov   + month_f, panel)
# Mooring days: log-linear OLS.
mA_moor_days <- lm_hac(log_mooring ~ gap_event_days + month_f, panel)
mA_moor_peak <- lm_hac(log_mooring ~ gap_peak_cov   + month_f, panel)
mA_moor_mean <- lm_hac(log_mooring ~ gap_mean_cov   + month_f, panel)

# -----------------------------------------------------------------------------
# 6. Model B - confound check: add year FE
# -----------------------------------------------------------------------------
# Both Bonaire's own event exposure and Bonaire tourism rise over 2022-2025
# (post-COVID recovery coinciding with growing regional sargassum). If the
# apparent effects are this shared trend, they should weaken or vanish once
# year fixed effects absorb the between-year movement.
# Yacht count: negative binomial (log link).
mB_cnt_days  <- nb_hac(yacht_n     ~ gap_event_days + month_f + year_f, panel)
mB_cnt_peak  <- nb_hac(yacht_n     ~ gap_peak_cov   + month_f + year_f, panel)
mB_cnt_mean  <- nb_hac(yacht_n     ~ gap_mean_cov   + month_f + year_f, panel)
# Mooring days: log-linear OLS.
mB_moor_days <- lm_hac(log_mooring ~ gap_event_days + month_f + year_f, panel)
mB_moor_peak <- lm_hac(log_mooring ~ gap_peak_cov   + month_f + year_f, panel)
mB_moor_mean <- lm_hac(log_mooring ~ gap_mean_cov   + month_f + year_f, panel)

# -----------------------------------------------------------------------------
# 7. Print HAC coefficient tables (sargassum terms only) to the console
# -----------------------------------------------------------------------------
report <- function(tag, fit) {
  keep <- grep("bon_|bar_|gap_", rownames(fit$ct))
  cat("---", tag, "---\n"); print(round(fit$ct[keep, , drop = FALSE], 4)); cat("\n")
}
cat("################  YACHT COUNT (negative binomial)  ################\n\n")
report("A  gap event_days (month FE)", mA_cnt_days)
report("A  gap peak_cov  (month FE)",  mA_cnt_peak)
report("A  gap mean_cov  (month FE)",  mA_cnt_mean)
report("B  gap event_days + year FE",  mB_cnt_days)
report("B  gap peak_cov  + year FE",   mB_cnt_peak)
report("B  gap mean_cov  + year FE",   mB_cnt_mean)
cat("################  MOORING DAYS (log-linear OLS)  ################\n\n")
report("A  gap event_days (month FE)", mA_moor_days)
report("A  gap peak_cov  (month FE)",  mA_moor_peak)
report("A  gap mean_cov  (month FE)",  mA_moor_mean)
report("B  gap event_days + year FE",  mB_moor_days)
report("B  gap peak_cov  + year FE",   mB_moor_peak)
report("B  gap mean_cov  + year FE",   mB_moor_mean)

# -----------------------------------------------------------------------------
# 8. Tidy comparison tables (one per outcome; HAC SEs)
# -----------------------------------------------------------------------------
library(flextable)   # modelsummary uses this for Word output (pulls in officer)

ms_ft <- function(mods, title) {
  modelsummary(
    mods,
    vcov     = lapply(mods, function(m) NeweyWest(m, lag = 3, prewhite = FALSE)),
    coef_map = gap_map,
    gof_omit = "AIC|BIC|Log.Lik|RMSE",
    stars    = c("*" = .05, "**" = .01, "***" = .001),
    title    = title,
    notes    = paste("Relative gap = Barbados minus Bonaire (detected_events.csv).",
                     "Month-of-year FE in all columns; B adds year FE.",
                     "Newey-West HAC SE, lag 3."),
    output   = "flextable"
  )
}

ft_count <- ms_ft(
  list("A: event days"       = mA_cnt_days$model,
       "A: peak cov"         = mA_cnt_peak$model,
       "A: mean cov"         = mA_cnt_mean$model,
       "B: event days +yrFE" = mB_cnt_days$model,
       "B: peak cov +yrFE"   = mB_cnt_peak$model,
       "B: mean cov +yrFE"   = mB_cnt_mean$model),
  "Table X. Displacement test - yacht count (negative binomial, HAC lag 3)"
)

ft_moor <- ms_ft(
  list("A: event days"       = mA_moor_days$model,
       "A: peak cov"         = mA_moor_peak$model,
       "A: mean cov"         = mA_moor_mean$model,
       "B: event days +yrFE" = mB_moor_days$model,
       "B: peak cov +yrFE"   = mB_moor_peak$model,
       "B: mean cov +yrFE"   = mB_moor_mean$model),
  "Table Y. Displacement test - log mooring days (OLS, HAC lag 3)"
)

# Both tables, each with its caption, in one .docx (written to getwd())
save_as_docx(
  `Yacht count (negative binomial)` = ft_count,
  `Mooring days (log-linear OLS)`   = ft_moor,
  path = "table_displacement.docx")