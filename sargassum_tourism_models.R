# =============================================================================
# Sargassum -> Tourism: distributed models for Bonaire
#   Outcomes : yacht mooring days, yacht count, cruise pax, aviation arrivals
#   Exposure : event-filtered intensity / days / peak coverage
#   Controls : month-of-year FE, GDP (84769NED) -- year FE as robustness
#   Sample   : monthly, COVID years (2020-2021) dropped
# =============================================================================


# .　 . • ☆ . ° .• °:. *₊ ° . ☆ SETUP .　 . • ☆ . ° .• °:. *₊ ° . ☆

library(dplyr)
library(lubridate)
library(ggplot2)
library(sandwich)   # NeweyWest()
library(lmtest)     # coeftest()
library(MASS)       # glm.nb()


# .　 . • ☆ . ° .• °:. *₊ ° . ☆ LOAD .　 . • ☆ . ° .• °:. *₊ ° . ☆

events <- read.csv("detected_events.csv")
events <- subset(events, island == "Bonaire")
events$start <- as.Date(events$start)
events$end   <- as.Date(events$end)
cat("Bonaire events:", nrow(events), "from",
    as.character(min(events$start)), "to",
    as.character(max(events$end)), "\n")

daily <- read.csv("sargassum_csvs/bonaire_coastal_all.csv")
daily$date <- as.Date(daily$date)


# .　 . • ☆ . ° .• °:. *₊ ° . ☆ EVENT-FILTERED DAILY EXPOSURE .　 . • ☆ . ° .• °:. *₊ ° . ☆

# Mark each day TRUE if it falls inside any detected event
daily$in_event <- FALSE
for (i in seq_len(nrow(events))) {
  hit <- daily$date >= events$start[i] & daily$date <= events$end[i]
  daily$in_event[hit] <- TRUE
}
cat("Daily rows inside an event:", sum(daily$in_event), "of", nrow(daily), "\n")

# Non-event days get zero exposure on every measure
daily$cov_event  <- ifelse(daily$in_event, daily$coverage_frac, 0)
daily$pos_event  <- daily$in_event & !is.na(daily$max_fai) & daily$max_fai > 0


# .　 . • ☆ . ° .• °:. *₊ ° . ☆ AGGREGATE TO MONTHLY .　 . • ☆ . ° .• °:. *₊ ° . ☆

daily$ym <- floor_date(daily$date, "month")

sarg_m <- daily |>
  group_by(ym) |>
  summarise(
    intensity = mean(cov_event, na.rm = TRUE),   # mean coverage_frac; non-event days are 0
    days_evt  = sum(pos_event,  na.rm = TRUE),   # event-positive days in the month
    peak_cov  = max(cov_event,  na.rm = TRUE),   # peak coverage in the month
    n_satobs  = sum(n_observed > 0, na.rm = TRUE)
  ) |>
  ungroup()

summary(sarg_m[, c("intensity", "days_evt", "peak_cov")])


# .　 . • ☆ . ° .• °:. *₊ ° . ☆ OUTCOMES .　 . • ☆ . ° .• °:. *₊ ° . ☆

yacht <- read.csv("cbs_raw/85015NED.csv")
yacht <- subset(yacht, Perioden_freq == "M")
yacht$ym <- as.Date(yacht$Perioden_Date)
yacht <- yacht[, c("ym", "Jachten_2", "GemiddeldAantalLigdagen_3")]
names(yacht) <- c("ym", "y_n", "y_moor")

cruise <- read.csv("cbs_raw/85007NED.csv")
cruise <- subset(cruise, Perioden_freq == "M")
cruise$ym <- as.Date(cruise$Perioden_Date)
cruise <- cruise[, c("ym", "Cruisepassagiers_1")]
names(cruise) <- c("ym", "cr_pax")

# Aviation: multi-airport dataset; verify Bonaire's label before filtering
aviation <- read.csv("cbs_raw/82332NED.csv")
aviation <- subset(aviation, Perioden_freq == "M")
cat("Airport labels in 82332NED:\n")
print(unique(aviation$LuchthavensCaribischNederland_label))
# ADJUST the regex below to match the actual Bonaire label
aviation <- aviation[grepl("Bonaire|Flamingo",
                           aviation$LuchthavensCaribischNederland_label,
                           ignore.case = TRUE), ]
aviation$ym <- as.Date(aviation$Perioden_Date)
aviation <- aviation[, c("ym", "AangekomenPassagiers_3")]
names(aviation) <- c("ym", "av_arr")


# .　 . • ☆ . ° .• °:. *₊ ° . ☆ MERGE & PAD .　 . • ☆ . ° .• °:. *₊ ° . ☆

# Use outcomes as the calendar backbone (they go back to 2012-14)
all_months <- seq(min(c(yacht$ym, cruise$ym, aviation$ym)),
                  max(c(yacht$ym, cruise$ym, aviation$ym)),
                  by = "month")
panel <- data.frame(ym = all_months)

panel <- merge(panel, sarg_m,   by = "ym", all.x = TRUE)
panel <- merge(panel, yacht,    by = "ym", all.x = TRUE)
panel <- merge(panel, cruise,   by = "ym", all.x = TRUE)
panel <- merge(panel, aviation, by = "ym", all.x = TRUE)

# Months before satellite series: by convention, exposure = 0 (no events).
# Inspect this carefully if any pre-2021 sargassum events are suspected.
panel$intensity[is.na(panel$intensity)] <- 0
panel$days_evt[is.na(panel$days_evt)]   <- 0
panel$peak_cov[is.na(panel$peak_cov)]   <- 0

panel$year  <- year(panel$ym)
panel$month <- month(panel$ym)
panel <- panel[order(panel$ym), ]




# .　 . • ☆ . ° .• °:. *₊ ° . ☆ DROP COVID & SET UP FACTORS .　 . • ☆ . ° .• °:. *₊ ° . ☆

panel_m <- subset(panel, (year %in% c(2022,2023,2024,2025,2026)))
panel_m$month_f <- factor(panel_m$month, levels = 1:12)
panel_m$year_f  <- factor(panel_m$year)
cat("Estimation sample rows:", nrow(panel_m), "\n")


# .　 . • ☆ . ° .• °:. *₊ ° . ☆ EDA .　 . • ☆ . ° .• °:. *₊ ° . ☆

p_expo <- ggplot(panel_m, aes(x = ym, y = intensity)) +
  geom_col(fill = "#E64B35", alpha = 0.8) +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  labs(title = "Monthly sargassum intensity (event-filtered)",
       x = NULL, y = "Mean coverage_frac") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_expo)

par(mfrow = c(2, 2))
plot(panel_m$ym, panel_m$y_moor, type = "l", main = "Yacht mooring days (avg)",
     xlab = NULL, ylab = "days")
plot(panel_m$ym, panel_m$y_n,    type = "l", main = "Yacht count",
     xlab = NULL, ylab = "yachts")
plot(panel_m$ym, panel_m$cr_pax, type = "l", main = "Cruise pax (x1000)",
     xlab = NULL, ylab = "x1000")
plot(panel_m$ym, panel_m$av_arr, type = "l", main = "Aviation arrivals",
     xlab = NULL, ylab = "passengers")
par(mfrow = c(1, 1))


panel_m %>%
  ggplot(aes(x = y_n)) + geom_histogram(bins = 30) + theme_bw() + 
  labs(x = "Monthly count of yacht", y = NULL)

panel_m %>%
  ggplot(aes(x = y_moor)) + geom_histogram(bins = 30) + theme_bw()+ 
  labs(x = "Monthly average of yacht mooring days", y = NULL)

panel_m %>%
  ggplot(aes(x = av_arr)) + geom_histogram(bins = 30) + theme_bw()+ 
  labs(x = "Monthly count of passengers arriving by flight", y = NULL)

panel_m %>%
  ggplot(aes(x = cr_pax)) + geom_histogram(bins = 30) + theme_bw()+ 
  labs(x = "Monthly count of passengers arriving by cruise ship", y = NULL)


library(AER)
rd <- glm(y_n ~ ., data = panel_m, family = poisson)
dispersiontest(rd,trafo=1)
rd <- glm(y_moor ~ ., data = panel_m, family = poisson)
dispersiontest(rd,trafo=1)
rd <- glm(av_arr ~ ., data = panel_m, family = poisson)
dispersiontest(rd,trafo=1)
rd <- glm(cr_pax ~ ., data = panel_m, family = poisson)
dispersiontest(rd,trafo=1)

panel_m$cr_pax<-panel_m$cr_pax*1000


# .　 . • ☆ . ° .• °:. *₊ ° . ☆ autocorrelation.　 . • ☆ . ° .• °:. *₊ ° . ☆

# ---- 1. Define the grid. EDIT to match your column names in panel_m ----
outcomes <- c(
  "Yacht mooring days" = "log(y_moor)",
  "Yacht count"        = "y_n",
  "Cruise passengers"  = "log(cr_pax+1)",
  "Aviation arrivals"  = "log(av_arr)"
)
exposures <- c(
  "Peak coverage"  = "peak_cov",
  "Days in event"  = "days_evt",
  "Mean intensity" = "intensity"
)

# ---- 2. Fit every model, collect residuals vs the exposure regressor ----
resid_long <- do.call(rbind, lapply(names(outcomes), function(oc) {
  do.call(rbind, lapply(names(exposures), function(ex) {
    f  <- reformulate(c(exposures[[ex]], "month_f", "year_f"),
                      response = outcomes[[oc]])
    m  <- lm(f, data = panel_m)
    mf <- model.frame(m)                 # rows actually used (NA-safe)
    data.frame(outcome = oc, exposure = ex,
               xval = mf[[exposures[[ex]]]], resi = resid(m))
  }))
}))
resid_long$outcome  <- factor(resid_long$outcome,  levels = names(outcomes))
resid_long$exposure <- factor(resid_long$exposure, levels = names(exposures))

# ---- 3. 12-panel residuals-vs-exposure grid ----
p_grid <- ggplot(resid_long, aes(xval, resi)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_point(colour = "blue", alpha = 0.55, size = 1.3) +
  facet_grid(outcome ~ exposure, scales = "free") +
  labs(x = "Exposure metric value", y = "OLS residual") +
  theme_bw(base_size = 10)

p_grid
ggsave("appendix_resid_grid.pdf", p_grid, width = 9, height = 11, units = "in")


# .　 . • ☆ . ° .• °:. *₊ ° . ☆ MAIN MODELS: Year FE control .　 . • ☆ . ° .• °:. *₊ ° . ☆

models <- list(
  "Mooring | intensity" = lm(log(y_moor)     ~ intensity + month_f + year_f, panel_m),
  "Mooring | days"      = lm(log(y_moor)     ~ days_evt  + month_f + year_f, panel_m),
  "Mooring | peak"      = lm(log(y_moor)     ~ peak_cov  + month_f + year_f, panel_m),
  
  "YachtN  | intensity" = glm.nb(y_n         ~ intensity + month_f + year_f, panel_m),
  "YachtN  | days"      = glm.nb(y_n         ~ days_evt  + month_f + year_f, panel_m),
  "YachtN  | peak"      = glm.nb(y_n         ~ peak_cov  + month_f + year_f, panel_m),
  
  "Cruise  | intensity" = lm(log(cr_pax + 1) ~ intensity + month_f + year_f, panel_m),
  "Cruise  | days"      = lm(log(cr_pax + 1) ~ days_evt  + month_f + year_f, panel_m),
  "Cruise  | peak"      = lm(log(cr_pax + 1) ~ peak_cov  + month_f + year_f, panel_m),
  
  "Avion   | intensity" = lm(log(av_arr)     ~ intensity + month_f + year_f, panel_m),
  "Avion   | days"      = lm(log(av_arr)     ~ days_evt  + month_f + year_f, panel_m),
  "Avion   | peak"      = lm(log(av_arr)     ~ peak_cov  + month_f + year_f, panel_m)
)

# HAC vcov for each model (Newey-West, lag 3)
vcov_list <- lapply(models, function(m) {
  NeweyWest(m, lag = 3, prewhite = FALSE)
})



# .　 . • ☆ . ° .• °:. *₊ ° . ☆ RESULTS TABLE .　 . • ☆ . ° .• °:. *₊ ° . ☆

library(modelsummary)

# Rename exposure terms to a common label so they stack on one row in the table
coef_map <- c(
  "intensity" = "Sargassum",
  "days_evt"  = "Sargassum",
  "peak_cov"  = "Sargassum"
)

# One table per outcome keeps each one readable (3 cols instead of 12).
outcomes <- c("Mooring", "YachtN", "Cruise", "Avion")
for (oc in outcomes) {
  idx <- grep(paste0("^", oc), names(models))
  modelsummary(
    models[idx],
    vcov = vcov_list[idx],
    coef_map = coef_map,
    coef_omit = "month_f|year_f|Intercept|theta",
    stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
    gof_omit = "AIC|BIC|RMSE|Log.Lik|F",
    notes = "HAC (Newey-West, lag 3) standard errors. Month-of-year and year FE included.",
    title = paste0("Sargassum effect on ", oc),
    output = paste0("table_", tolower(oc), ".docx")
  )
}

# .　 . • ☆ . ° .• °:. *₊ ° . ☆ COEF PLOT WITH SIGNIFICANCE .　 . • ☆ . ° .• °:. *₊ ° . ☆

# Add p-values to results_df (rebuild from cell_df above to keep them aligned)
results_df <- merge(
  results_df,
  cell_df[, c("outcome", "exposure", "p")],
  by = c("outcome", "exposure")
)


results_df$outcome <- factor(results_df$outcome,
                             levels = c("Mooring", "YachtN", "Avion", "Cruise"))
results_df$exposure <- factor(results_df$exposure,
                              levels = c("intensity", "days", "peak"))

ggplot(results_df, aes(x = estimate, y = outcome,
                       color = exposure, shape = sig, fill = exposure)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                 position = position_dodge(width = 0.55),
                 linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.55), size = 3.2, stroke = 1) +
  geom_text(aes(label = siglab),
            position = position_dodge(width = 0.55),
            vjust = -1.0, size = 4.5, show.legend = FALSE) +
  scale_color_manual(values = c(intensity = "#E64B35",
                                days      = "#3C5488",
                                peak      = "#00A087")) +
  scale_fill_manual(values  = c(intensity = "#E64B35",
                                days      = "#3C5488",
                                peak      = "#00A087")) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 21),
                     guide  = "none") +    # hide shape legend; fill carries it
  # solid fill when significant, white fill when not
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0)) +
  labs(title = "Sargassum effect on Bonaire tourism outcomes",
       subtitle = "1 SD increase in exposure; 95% HAC CI. Solid = p<0.05, hollow = n.s. Stars: * <0.1, ** <0.05, *** <0.01",
       x = "Effect on log(outcome)", y = NULL,
       color = "Exposure", fill = "Exposure") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "bottom")






# .　 . • ☆ . ° .• °:. *₊ ° . ☆ Robust tests. .　 . • ☆ . ° .• °:. *₊ ° . ☆


vcov_list2 <- lapply(models, function(m) {
  NeweyWest(m, lag = 2, prewhite = FALSE)
})

vcov_list4 <- lapply(models, function(m) {
  NeweyWest(m, lag = 4, prewhite = FALSE)
})

vcov_list6 <- lapply(models, function(m) {
  NeweyWest(m, lag = 6, prewhite = FALSE)
})



# .　 . • ☆ . ° .• °:. *₊ ° . ☆ RESULTS TABLE .　 . • ☆ . ° .• °:. *₊ ° . ☆

models %>% modelsummary(
  vcov = vcov_list6,
  coef_map = coef_map,
  coef_omit = "month_f|year_f|Intercept|theta",
  stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  gof_omit = "AIC|BIC|RMSE|Log.Lik|F",
  notes = "HAC (Newey-West, lag 3) standard errors. Month-of-year and year FE included.",
  title = paste0("Sargassum effect on everything")
)

# robust test passed



# .　 . • ☆ . ° .• °:. *₊ ° . ☆ PARTIAL RESIDUAL CROSS-PLOTS .　 . • ☆ . ° .• °:. *₊ ° . ☆

library(ggplot2)
library(ggrepel)   # install.packages("ggrepel") if needed

# Residualize a series against month FE and year FE on the estimation sample.
# Returns a vector aligned to panel_m row order, NA where the input is NA.
resid_fe <- function(y, data = panel_m) {
  ok <- complete.cases(y, data$month_f, data$year_f)
  out <- rep(NA_real_, length(y))
  out[ok] <- residuals(lm(y[ok] ~ data$month_f[ok] + data$year_f[ok]))
  out
}

# Build a plotting frame with residualized exposure and residualized outcomes.
# For yacht count (NB outcome), use log(y_n) so residuals are on a comparable
# multiplicative scale; doesn't have to match the GLM exactly because the
# purpose here is diagnostic, not estimation.
cross_df <- data.frame(
  ym         = panel_m$ym,
  year       = panel_m$year,
  month      = panel_m$month,
  # Residualized exposures (one per measure)
  r_intensity = resid_fe(panel_m$intensity),
  r_days      = resid_fe(panel_m$days_evt),
  r_peak      = resid_fe(panel_m$peak_cov),
  # Residualized outcomes
  r_moor   = resid_fe(log(panel_m$y_moor)),
  r_yn     = resid_fe(log(panel_m$y_n)),
  r_cruise = resid_fe(log(panel_m$cr_pax + 1)),
  r_avion  = resid_fe(log(panel_m$av_arr))
)

# Long format for faceting: one row per (outcome x exposure) cell
library(tidyr)
long <- cross_df |>
  pivot_longer(cols = starts_with("r_moor"):starts_with("r_avion"),
               names_to = "outcome", values_to = "y_resid") |>
  pivot_longer(cols = c(r_intensity, r_days, r_peak),
               names_to = "exposure", values_to = "x_resid")

# Pretty labels
long$outcome  <- factor(long$outcome,
                        levels = c("r_moor", "r_yn", "r_cruise", "r_avion"),
                        labels = c("Yacht mooring days", "Yacht count",
                                   "Cruise passengers",  "Aviation arrivals"))
long$exposure <- factor(long$exposure,
                        levels = c("r_intensity", "r_days", "r_peak"),
                        labels = c("Intensity", "Days in event", "Peak coverage"))

# Label only points with extreme leverage (|x_resid| in top decile)
long$ymlabel <- ifelse(
  !is.na(long$x_resid) &
    abs(long$x_resid) >= quantile(abs(long$x_resid), 0.9, na.rm = TRUE),
  format(long$ym, "%Y-%m"),
  NA_character_
)

ggplot(long, aes(x = x_resid, y = y_resid)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_point(aes(color = factor(year)), size = 2, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  geom_text_repel(aes(label = ymlabel), size = 3,
                  max.overlaps = Inf, segment.alpha = 0.4, na.rm = TRUE) +
  facet_grid(outcome ~ exposure, scales = "free") +
  scale_color_brewer(palette = "Set2", name = "Year") +
  labs(title = "Partial-residual cross-plots: sargassum vs. tourism outcomes",
       subtitle = "Both axes residualized against month-of-year and year fixed effects. Slope = model coefficient.",
       x = "Sargassum exposure (residualized)",
       y = "Outcome, log scale (residualized)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))