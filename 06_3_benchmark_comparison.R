# ============================================================
# Phase-stratified error comparison: proposed model vs SARIMA
# Independent holdout period, weeks 169-200
#run this code after file 06_2_bechmaek_independent.R
# ============================================================

library(dplyr)

decline_df <- data.frame(
  week     = test_full$time_index,
  observed = obs_test,
  proposed = pred_proposed,
  sarima   = pred_sarima,
  ar_only  = pred_ar_only
) %>%
  arrange(week)

# Previous observed cases:
# Week 169 uses the actual observed value from week 168;
# subsequent weeks use the preceding holdout observation.
previous_cases <- c(
  df$Dengue_cases[df$time_index == min(decline_df$week) - 1],
  head(decline_df$observed, -1)
)

decline_df <- decline_df %>%
  mutate(
    previous_observed = previous_cases,
    change = observed - previous_observed,
    
    phase = case_when(
      change < 0 ~ "Decreasing",
      change > 0 ~ "Increasing",
      change == 0 ~ "Stable",
      TRUE ~ NA_character_
    ),
    
    # Signed error: observed - predicted
    # Negative values indicate overprediction
    error_proposed = observed - proposed,
    error_sarima   = observed - sarima,
    error_ar       = observed - ar_only,
    
    AE_proposed = abs(error_proposed),
    AE_sarima   = abs(error_sarima),
    AE_ar       = abs(error_ar),
    
    # Positive values indicate smaller SARIMA absolute error
    SARIMA_advantage = AE_proposed - AE_sarima
  )


# ------------------------------------------------------------
# Phase-specific performance
# ------------------------------------------------------------

phase_summary <- decline_df %>%
  filter(!is.na(phase)) %>%
  group_by(phase) %>%
  summarise(
    N = n(),
    Proposed_MAE = mean(AE_proposed, na.rm = TRUE),
    SARIMA_MAE   = mean(AE_sarima, na.rm = TRUE),
    AR_only_MAE  = mean(AE_ar, na.rm = TRUE),
    Proposed_mean_error = mean(error_proposed, na.rm = TRUE),
    SARIMA_mean_error   = mean(error_sarima, na.rm = TRUE),
    SARIMA_better_n = sum(AE_sarima < AE_proposed, na.rm = TRUE),
    SARIMA_better_pct =
      100 * mean(AE_sarima < AE_proposed, na.rm = TRUE),
    .groups = "drop"
  )

print(phase_summary)


# ------------------------------------------------------------
# Declining-week overprediction
# ------------------------------------------------------------

declining_weeks <- decline_df %>%
  filter(phase == "Decreasing")

decline_bias <- declining_weeks %>%
  summarise(
    N = n(),
    
    Proposed_mean_error =
      mean(error_proposed, na.rm = TRUE),
    
    Proposed_median_error =
      median(error_proposed, na.rm = TRUE),
    
    SARIMA_mean_error =
      mean(error_sarima, na.rm = TRUE),
    
    SARIMA_median_error =
      median(error_sarima, na.rm = TRUE),
    
    Proposed_overprediction_n =
      sum(error_proposed < 0, na.rm = TRUE),
    
    Proposed_overprediction_pct =
      100 * mean(error_proposed < 0, na.rm = TRUE),
    
    SARIMA_overprediction_n =
      sum(error_sarima < 0, na.rm = TRUE),
    
    SARIMA_overprediction_pct =
      100 * mean(error_sarima < 0, na.rm = TRUE)
  )

print(decline_bias)


# ------------------------------------------------------------
# Paired comparison during declining weeks
# ------------------------------------------------------------

wilcox_decline <- wilcox.test(
  declining_weeks$AE_proposed,
  declining_weeks$AE_sarima,
  paired = TRUE,
  exact = FALSE
)

print(wilcox_decline)


# ------------------------------------------------------------
# SARIMA advantage by epidemic phase
# ------------------------------------------------------------

advantage_by_phase <- decline_df %>%
  filter(phase %in% c("Increasing", "Decreasing")) %>%
  group_by(phase) %>%
  summarise(
    N = n(),
    Mean_advantage =
      mean(SARIMA_advantage, na.rm = TRUE),
    Median_advantage =
      median(SARIMA_advantage, na.rm = TRUE),
    SARIMA_better_pct =
      100 * mean(SARIMA_advantage > 0, na.rm = TRUE),
    .groups = "drop"
  )

print(advantage_by_phase)