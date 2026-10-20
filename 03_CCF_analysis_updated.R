# ============================================================
# RAW CROSS-CORRELATION FUNCTION (CCF)
# DENGUE FORECASTING IN BANGLADESH
#
# Training: Weeks 1-168
# Holdout:  Weeks 169-200 (EXCLUDED)
#
# Meteorological variables:
#   1. Humidity
#   2. Precipitation
#   3. Temperature
#
# Lag range: 0-24 weeks
# Blue dashed lines: Approximate 95% reference bounds
# ============================================================

library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load and order data
# ------------------------------------------------------------

main_data_path <- paste0(
  "/Users/akidulhaque/Downloads/",
  "MyBDweeklyFinaldata.csv"
)

data1 <- read.csv(
  main_data_path,
  stringsAsFactors = FALSE
) %>%
  arrange(Week)

# ------------------------------------------------------------
# 2. Select training data only
# ------------------------------------------------------------

train <- data1 %>%
  filter(Week >= 1 & Week <= 168)

holdout <- data1 %>%
  filter(Week >= 169 & Week <= 200)

# Verify chronological split
stopifnot(
  nrow(data1) == 200,
  nrow(train) == 168,
  nrow(holdout) == 32,
  identical(as.numeric(train$Week), as.numeric(1:168)),
  identical(as.numeric(holdout$Week), as.numeric(169:200))
)

required_cols <- c(
  "Dengue_cases",
  "Temperature",
  "Precipitation",
  "Humidity"
)

stopifnot(
  all(required_cols %in% names(train)),
  !anyNA(train[, required_cols])
)

cat("Training observations:", nrow(train), "\n")
cat("Training weeks:", min(train$Week),
    "-", max(train$Week), "\n")
cat("Holdout excluded: TRUE\n")

# ------------------------------------------------------------
# 3. Set CCF parameters
# ------------------------------------------------------------

max_lag <- 24
peak_fraction <- 0.85

# Approximate 95% reference bounds
n_train <- nrow(train)

ccf_bound <- 1.96 / sqrt(n_train)

cat(
  "Approximate 95% CCF bounds: +/-",
  round(ccf_bound, 1), "\n"
)

# ------------------------------------------------------------
# 4. Function to calculate raw CCF
# ------------------------------------------------------------
# Positive reported lag k means:
# Meteorological variable at t-k is associated
# with dengue cases at week t.
#
# R's ccf() uses the opposite lag sign for
# ccf(exposure, outcome), so we reverse it.
# ------------------------------------------------------------

run_raw_ccf <- function(
    exposure,
    outcome,
    max_lag = 24
) {
  
  cc <- ccf(
    x = exposure,
    y = outcome,
    lag.max = max_lag,
    plot = FALSE,
    na.action = na.fail
  )
  
  result <- data.frame(
    Lag = -as.numeric(cc$lag),
    Correlation = as.numeric(cc$acf)
  )
  
  result %>%
    filter(
      Lag >= 0,
      Lag <= max_lag
    ) %>%
    arrange(Lag)
}

# ------------------------------------------------------------
# 5. Calculate CCF for each variable
# ------------------------------------------------------------

ccf_humid <- run_raw_ccf(
  train$Humidity,
  train$Dengue_cases,
  max_lag
)

ccf_precip <- run_raw_ccf(
  train$Precipitation,
  train$Dengue_cases,
  max_lag
)

ccf_temp <- run_raw_ccf(
  train$Temperature,
  train$Dengue_cases,
  max_lag
)

# ------------------------------------------------------------
# 6. Combine CCF results
# ------------------------------------------------------------

ccf_all <- bind_rows(
  
  ccf_humid %>%
    mutate(Variable = "Humidity"),
  
  ccf_precip %>%
    mutate(Variable = "Precipitation"),
  
  ccf_temp %>%
    mutate(Variable = "Temperature")
  
) %>%
  mutate(
    Variable = factor(
      Variable,
      levels = c(
        "Humidity",
        "Precipitation",
        "Temperature"
      )
    )
  )

# ------------------------------------------------------------
# 7. Candidate lag windows (90% of peak)
# ------------------------------------------------------------

get_candidate_window <- function(
    df,
    prop_of_peak = 0.90
) {
  
  df <- df %>%
    filter(Lag >= 1) %>%
    arrange(Lag)
  
  peak_index <- which.max(df$Correlation)
  
  peak_lag <- df$Lag[peak_index]
  peak_corr <- df$Correlation[peak_index]
  
  if (peak_corr <= 0) {
    return(data.frame(
      Peak_lag = peak_lag,
      Peak_correlation = peak_corr,
      Window_start = NA_integer_,
      Window_end = NA_integer_
    ))
  }
  
  threshold <- prop_of_peak * peak_corr
  
  candidate_lags <- df$Lag[
    df$Correlation >= threshold
  ]
  
  # Find contiguous sequences
  runs <- split(
    candidate_lags,
    cumsum(
      c(1, diff(candidate_lags) != 1)
    )
  )
  
  # Select the contiguous run containing the peak
  selected_run <- Filter(
    function(x) peak_lag %in% x,
    runs
  )[[1]]
  
  data.frame(
    Peak_lag = peak_lag,
    Peak_correlation = peak_corr,
    Window_start = min(selected_run),
    Window_end = max(selected_run)
  )
}

ccf_summary <- ccf_all %>%
  group_by(Variable) %>%
  group_modify(
    ~ get_candidate_window(
      .x,
      prop_of_peak = peak_fraction
    )
  ) %>%
  ungroup()

cat("\n====================================\n")
cat("TRAINING-ONLY CCF RESULTS\n")
cat("====================================\n")

print(ccf_summary, n = Inf)

# ------------------------------------------------------------
# 8. Publication theme
# ------------------------------------------------------------

theme_pub <- function(base_size = 13) {
  
  theme_bw(base_size = base_size) +
    
    theme(
      
      panel.grid = element_blank(),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.8
      ),
      
      strip.background = element_rect(
        fill = "grey90",
        colour = "black",
        linewidth = 0.5
      ),
      
      strip.text = element_text(
        face = "bold",
        colour = "black"
      ),
      
      axis.title = element_text(
        face = "bold"
      ),
      
      axis.text = element_text(
        colour = "black"
      )
    )
}

# ------------------------------------------------------------
# 9. Plot raw CCF with blue dashed lines
# ------------------------------------------------------------

p_ccf <- ggplot(
  ccf_all,
  aes(
    x = Lag,
    y = Correlation
  )
) +
  
  # Zero reference line
  geom_hline(
    yintercept = 0,
    colour = "black",
    linewidth = 0.5
  ) +
  
  # Upper and lower 95% reference bounds
  geom_hline(
    yintercept = c(-ccf_bound, ccf_bound),
    colour = "#2C7FB8",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  # Vertical CCF bars
  geom_segment(
    aes(
      xend = Lag,
      yend = 0
    ),
    colour = "black",
    linewidth = 0.7
  ) +
  
  # Points at the end of bars
  geom_point(
    colour = "black",
    size = 1.8
  ) +
  
  # Three panels in one row
  facet_wrap(
    ~ Variable,
    nrow = 1,
    scales = "fixed"
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 24, by = 2),
    limits = c(0, 24)
  ) +
  
  scale_y_continuous(
    breaks = seq(-0.6, 0.6, by = 0.2),
    labels = function(x) sprintf("%.1f", x),
    limits = c(-0.6, 0.6)
  )+
  
  labs(
    x = "Lag (weeks)",
    y = "Cross-correlation coefficient"
  ) +
  
  theme_pub(base_size = 13)

print(p_ccf)

# ------------------------------------------------------------
# 10. Save figure
# ------------------------------------------------------------

ggsave(
  filename = "CCF_training.png",
  plot = p_ccf,
  width = 10,
  height = 4,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "CCF_training.pdf",
  plot = p_ccf,
  width = 10,
  height = 4,
  units = "in"
)

# ------------------------------------------------------------
# 11. Export results
# ------------------------------------------------------------

write.csv(
  ccf_all,
  "CCF_training_correlations.csv",
  row.names = FALSE
)

write.csv(
  ccf_summary,
  "CCF_training_lag_windows.csv",
  row.names = FALSE
)

cat("\nCCF analysis completed successfully.\n")
cat("All calculations used weeks 1-168 only.\n")
cat("Holdout weeks 169-200 were excluded.\n")