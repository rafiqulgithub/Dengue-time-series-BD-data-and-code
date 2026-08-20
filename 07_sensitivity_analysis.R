
library(dplyr)
library(splines)
library(knitr)

df <-read.csv("MyBDweeklyFinaldata.csv",stringsAsFactors = FALSE)

df <- df %>%
  arrange(Week) %>%
  mutate(
    time_index = row_number(),
    trend = time_index,
    sin1 = sin(2 * pi * time_index / 52),
    cos1 = cos(2 * pi * time_index / 52)
  )

final_ar_order   <- 8
final_temp_lag   <- 9
final_precip_lag <- 7
final_humid_lag  <- 8

rmse <- function(o, p) sqrt(mean((o - p)^2, na.rm = TRUE))
mae  <- function(o, p) mean(abs(o - p), na.rm = TRUE)

make_model_data <- function(data, ar_order, temp_lag, precip_lag, humid_lag) {
  out <- data
  for (i in seq_len(ar_order)) {
    out[[paste0("lag_dengue_", i)]] <- log(dplyr::lag(out$Dengue_cases, i) + 1)
  }
  out$Temp_lag   <- dplyr::lag(out$Temperature, temp_lag)
  out$Precip_lag <- dplyr::lag(out$Precipitation, precip_lag)
  out$Humd_lag   <- dplyr::lag(out$Humidity, humid_lag)
  out
}

# knot_probs is now a parameter, defaulting to your baseline (25/50/75)
build_ns_basis <- function(train_x, knot_probs = c(0.25, 0.50, 0.75)) {
  train_x <- train_x[is.finite(train_x)]
  kn <- unique(as.numeric(stats::quantile(train_x, probs = knot_probs, na.rm = TRUE)))
  bk <- range(train_x, na.rm = TRUE)
  basis_train <- ns(train_x, knots = kn, Boundary.knots = bk)
  list(train_basis = basis_train, predict_fn = function(newx) predict(basis_train, newx))
}

add_spline_terms <- function(train_model, valid_model, knot_probs = c(0.25, 0.50, 0.75)) {
  temp_basis   <- build_ns_basis(train_model$Temp_lag,   knot_probs)
  precip_basis <- build_ns_basis(train_model$Precip_lag, knot_probs)
  humid_basis  <- build_ns_basis(train_model$Humd_lag,   knot_probs)
  
  train_temp <- as.data.frame(temp_basis$train_basis)
  train_prec <- as.data.frame(precip_basis$train_basis)
  train_humd <- as.data.frame(humid_basis$train_basis)
  
  colnames(train_temp) <- paste0("ns_temp_", seq_len(ncol(train_temp)))
  colnames(train_prec) <- paste0("ns_precip_", seq_len(ncol(train_prec)))
  colnames(train_humd) <- paste0("ns_humid_", seq_len(ncol(train_humd)))
  
  valid_temp <- as.data.frame(temp_basis$predict_fn(valid_model$Temp_lag))
  valid_prec <- as.data.frame(precip_basis$predict_fn(valid_model$Precip_lag))
  valid_humd <- as.data.frame(humid_basis$predict_fn(valid_model$Humd_lag))
  
  colnames(valid_temp) <- colnames(train_temp)
  colnames(valid_prec) <- colnames(train_prec)
  colnames(valid_humd) <- colnames(train_humd)
  
  train_full <- bind_cols(train_model, train_temp, train_prec, train_humd)
  valid_full <- bind_cols(valid_model, valid_temp, valid_prec, valid_humd)
  
  spline_terms <- c(colnames(train_temp), colnames(train_prec), colnames(train_humd))
  
  list(train_full = train_full, valid_full = valid_full, spline_terms = spline_terms)
}

# -----------------------------
# Generic fit + evaluate function, parameterized by split and knots
# -----------------------------
run_sensitivity_config <- function(train_end_week, test_start_week, knot_probs, label) {
  
  max_lag <- max(final_ar_order, final_temp_lag, final_precip_lag, final_humid_lag)
  
  train_raw <- df %>% filter(time_index <= train_end_week)
  test_raw  <- df %>% filter(time_index >= (train_end_week - max_lag + 1))
  
  train_lagged <- make_model_data(train_raw, final_ar_order, final_temp_lag, final_precip_lag, final_humid_lag)
  test_lagged  <- make_model_data(test_raw,  final_ar_order, final_temp_lag, final_precip_lag, final_humid_lag)
  
  ar_terms <- paste0("lag_dengue_", seq_len(final_ar_order))
  needed_cols <- c("Dengue_cases", ar_terms, "Temp_lag", "Precip_lag", "Humd_lag",
                   "trend", "sin1", "cos1")
  
  train_model <- train_lagged %>% filter(if_all(all_of(needed_cols), ~ !is.na(.)))
  test_model <- test_lagged %>%
    filter(if_all(all_of(needed_cols), ~ !is.na(.))) %>%
    filter(time_index >= test_start_week)
  
  spline_obj <- add_spline_terms(train_model, test_model, knot_probs)
  train_full <- spline_obj$train_full
  test_full  <- spline_obj$valid_full
  spline_terms <- spline_obj$spline_terms
  
  model_formula <- as.formula(
    paste("Dengue_cases ~", paste(ar_terms, collapse = " + "), "+",
          paste(spline_terms, collapse = " + "), "+", "trend + sin1 + cos1")
  )
  
  fit <- glm(model_formula, family = quasipoisson(link = "log"), data = train_full)
  
  train_pred <- predict(fit, newdata = train_full, type = "response")
  test_pred  <- predict(fit, newdata = test_full,  type = "response")
  train_sd   <- sd(train_full$Dengue_cases, na.rm = TRUE)
  
  data.frame(
    Config = label,
    N_train = nrow(train_full), N_test = nrow(test_full),
    Train_RMSE = rmse(train_full$Dengue_cases, train_pred),
    Test_RMSE  = rmse(test_full$Dengue_cases,  test_pred),
    Test_MAE   = mae(test_full$Dengue_cases,   test_pred),
    Test_SRMSE = rmse(test_full$Dengue_cases,  test_pred) / train_sd,
    Test_SMAE  = mae(test_full$Dengue_cases,   test_pred) / train_sd
  )
}

# ============================================================
# ARM 1: Split ratio sensitivity (84/16 baseline vs 80/20, 70/30)
# Knots held at baseline (25/50/75)
# ============================================================
split_results <- bind_rows(
  run_sensitivity_config(168, 169, c(0.25, 0.50, 0.75), "84/16 (baseline)"),
  run_sensitivity_config(160, 161, c(0.25, 0.50, 0.75), "80/20"),
  run_sensitivity_config(140, 141, c(0.25, 0.50, 0.75), "70/30")
)

cat("\n--- Sensitivity: train/test split ratio ---\n")
print(split_results)
kable(split_results, digits = 4, caption = "Sensitivity of model performance to train/test split ratio")

# ============================================================
# ARM 2: Knot placement sensitivity
# Baseline (25/50/75) vs 10/50/90 vs 20/50/80
# Split held at baseline (84/16)
# ============================================================
knot_results <- bind_rows(
  run_sensitivity_config(168, 169, c(0.25, 0.50, 0.75), "Knots: 25/50/75 (baseline)"),
  run_sensitivity_config(168, 169, c(0.10, 0.50, 0.90), "Knots: 10/50/90"),
  run_sensitivity_config(168, 169, c(0.20, 0.50, 0.80), "Knots: 20/50/80")
)

cat("\n--- Sensitivity: spline knot placement ---\n")
print(knot_results)
kable(knot_results, digits = 4, caption = "Sensitivity of model performance to spline knot placement")

# ============================================================
# ARM 3: Model family (quasi-Poisson vs negative binomial)
# NOT re-run here -- reused from the baseline benchmark comparison
# (evaluate_one_split() in the benchmark script), both fit at the
# 84/16 split with 25/50/75 knots:
#   Proposed quasi-Poisson: SRMSE = 0.0998
#   Negative binomial:      SRMSE = 0.1007
# ============================================================
cat("\n--- Sensitivity: model family (reused from baseline benchmark) ---\n")
family_results <- data.frame(
  Config = c("Quasi-Poisson (baseline)", "Negative binomial"),
  Test_SRMSE = c(0.0998, 0.1007)
)
print(family_results)
kable(family_results, digits = 4, caption = "Sensitivity of model performance to model family (reused from baseline benchmark comparison)")

# -----------------------------
# Save all results
# -----------------------------
write.csv(split_results, "sensitivity_split_ratio.csv", row.names = FALSE)
write.csv(knot_results,  "sensitivity_knot_placement.csv", row.names = FALSE)
write.csv(family_results, "sensitivity_model_family.csv", row.names = FALSE)

sd(df$Dengue_cases[169:200])
sd(df$Dengue_cases[161:200])
sd(df$Dengue_cases[141:200])