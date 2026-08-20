

# ==============================================================================

library(dplyr)
library(splines)
library(MASS)
library(forecast)
library(knitr)

# -----------------------------
# Output directory -- change to your local folder
# -----------------------------
out_dir <- "~/Downloads/benchmark_outputs_cv"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load data and create predictors
# ------------------------------------------------------------------------------
df <- read.csv("C:\\Users\\Rafiqul islam\\Desktop\\Reaearch\\Akib\\Dengue\\MyBDweeklyFinaldata.csv",
               stringsAsFactors = FALSE)

df <- df %>%
  arrange(Week) %>%
  mutate(
    time_index = row_number(),
    trend = time_index,
    sin1 = sin(2 * pi * time_index / 52),
    cos1 = cos(2 * pi * time_index / 52),
    Temp_lag   = dplyr::lag(Temperature, 9),
    Precip_lag = dplyr::lag(Precipitation, 7),
    Humd_lag   = dplyr::lag(Humidity, 8)
  )

for (i in 1:8) {
  df[[paste0("lag_dengue_", i)]] <- log(dplyr::lag(df$Dengue_cases, i) + 1)
}

ar_terms <- paste0("lag_dengue_", 1:8, collapse = " + ")

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae  <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)

calc_metrics <- function(obs, pred, train_sd) {
  data.frame(
    RMSE  = round(rmse(obs, pred), 2),
    MAE   = round(mae(obs, pred), 2),
    SRMSE = round(rmse(obs, pred) / train_sd, 3),
    SMAE  = round(mae(obs, pred) / train_sd, 3)
  )
}

build_ns_basis <- function(train_x) {
  x_clean <- train_x[is.finite(train_x)]
  knots <- as.numeric(quantile(x_clean, probs = c(0.25, 0.50, 0.75), na.rm = TRUE))
  knots <- unique(knots)
  boundary_knots <- range(x_clean, na.rm = TRUE)
  basis_train <- ns(train_x, knots = knots, Boundary.knots = boundary_knots)
  list(train_basis = basis_train,
       predict_fn  = function(new_x) predict(basis_train, new_x))
}

# ==============================================================================
# 2. Reusable function: fit all five models on one train/test split
# ==============================================================================
evaluate_one_split <- function(train_start, train_end, test_start, test_end,
                               split_label = NA) {
  
  train_data <- df %>%
    filter(time_index >= train_start, time_index <= train_end) %>%
    filter(
      !is.na(lag_dengue_8), !is.na(Temp_lag),
      !is.na(Precip_lag),   !is.na(Humd_lag)
    )
  
  test_data <- df %>%
    filter(time_index >= test_start, time_index <= test_end) %>%
    filter(
      !is.na(lag_dengue_8), !is.na(Temp_lag),
      !is.na(Precip_lag),   !is.na(Humd_lag)
    )
  
  if (nrow(train_data) < 20 || nrow(test_data) < 1) {
    warning("Split '", split_label, "' has insufficient rows; skipping.")
    return(NULL)
  }
  
  # Full development set for univariate SARIMA (no lag-based NA
  # dropping needed, since SARIMA only uses the case-count series)
  sarima_train_data <- df %>%
    filter(time_index >= train_start, time_index <= train_end)
  
  train_sd <- sd(train_data$Dengue_cases, na.rm = TRUE)
  test_sd<- sd(test_data$Dengue_cases, na.rm = TRUE)
  obs_test <- test_data$Dengue_cases
  
  # -----------------------------
  # Spline bases (fit on training data only)
  # -----------------------------
  temp_basis   <- build_ns_basis(train_data$Temp_lag)
  precip_basis <- build_ns_basis(train_data$Precip_lag)
  humid_basis  <- build_ns_basis(train_data$Humd_lag)
  
  train_splines <- data.frame(
    ns_temp   = as.matrix(temp_basis$train_basis),
    ns_precip = as.matrix(precip_basis$train_basis),
    ns_humid  = as.matrix(humid_basis$train_basis)
  )
  test_splines <- data.frame(
    ns_temp   = as.matrix(temp_basis$predict_fn(test_data$Temp_lag)),
    ns_precip = as.matrix(precip_basis$predict_fn(test_data$Precip_lag)),
    ns_humid  = as.matrix(humid_basis$predict_fn(test_data$Humd_lag))
  )
  colnames(train_splines) <- paste0("ns_col_", seq_len(ncol(train_splines)))
  colnames(test_splines)  <- colnames(train_splines)
  
  train_full <- bind_cols(train_data, train_splines)
  test_full  <- bind_cols(test_data, test_splines)
  
  spline_terms <- paste(colnames(train_splines), collapse = " + ")
  formula_proposed <- as.formula(paste("Dengue_cases ~", ar_terms, "+", spline_terms, "+ trend + sin1 + cos1"))
  formula_ar_only  <- as.formula(paste("Dengue_cases ~", ar_terms, "+ trend + sin1 + cos1"))
  
  # -----------------------------
  # Model 1: Proposed quasi-Poisson
  # -----------------------------
  fit_proposed  <- glm(formula_proposed, family = quasipoisson(link = "log"), data = train_full)
  pred_proposed <- predict(fit_proposed, newdata = test_full, type = "response")
  
  # -----------------------------
  # Model 2: Negative binomial
  # -----------------------------
  fit_nb <- tryCatch(
    MASS::glm.nb(formula_proposed, data = train_full),
    error = function(e) {
      warning("Negative binomial failed for ", split_label, ": ", conditionMessage(e))
      NULL
    }
  )
  pred_nb <- if (!is.null(fit_nb)) {
    predict(fit_nb, newdata = test_full, type = "response")
  } else {
    rep(NA_real_, nrow(test_full))
  }
  
  # -----------------------------
  # Model 3: AR(8) quasi-Poisson (no meteorology)
  # -----------------------------
  fit_ar_only  <- glm(formula_ar_only, family = quasipoisson(link = "log"), data = train_full)
  pred_ar_only <- predict(fit_ar_only, newdata = test_full, type = "response")
  
  # -----------------------------
  # Model 4: Sequential one-step-ahead SARIMA
  # Structure/coefficients fixed from the FULL development period for
  # this fold, then sequentially reapplied to the growing observed
  # history (true values revealed after each one-step forecast).
  # -----------------------------
  ts_train <- ts(sarima_train_data$Dengue_cases, frequency = 52)
  fit_sarima_train <- tryCatch(
    forecast::auto.arima(ts_train, seasonal = TRUE, stepwise = FALSE, approximation = FALSE),
    error = function(e) {
      warning("SARIMA structure selection failed for ", split_label, ": ", conditionMessage(e))
      NULL
    }
  )
  
  pred_sarima <- rep(NA_real_, nrow(test_data))
  
  if (!is.null(fit_sarima_train)) {
    history_vec <- sarima_train_data$Dengue_cases
    
    for (i in seq_len(nrow(test_data))) {
      history_ts <- ts(history_vec, frequency = 52)
      
      fit_iter <- tryCatch(
        forecast::Arima(history_ts, model = fit_sarima_train),
        error = function(e) {
          warning("SARIMA reapplication failed for ", split_label, " at step ", i,
                  ": ", conditionMessage(e))
          NULL
        }
      )
      
      pred_sarima[i] <- if (!is.null(fit_iter)) {
        as.numeric(forecast::forecast(fit_iter, h = 1)$mean[1])
      } else {
        NA_real_
      }
      
      history_vec <- c(history_vec, test_data$Dengue_cases[i])
    }
  }
  
  # -----------------------------
  # Model 5: Seasonal naive (Y_t = Y_{t-52})
  # -----------------------------
  pred_snaive <- sapply(test_full$time_index, function(ti) {
    val <- df$Dengue_cases[df$time_index == (ti - 52)]
    if (length(val) == 1) val else NA_real_
  })
  
  # -----------------------------
  # Ensure all five models are evaluated on identical weeks.
  # If any model fails to produce a prediction for a given week
  # (e.g. SARIMA reapplication fails, NB fails to converge), that
  # week is dropped for ALL models so the comparison stays fair.
  # -----------------------------
  prediction_df <- data.frame(
    observed           = obs_test,
    proposed           = pred_proposed,
    negative_binomial  = pred_nb,
    ar_only            = pred_ar_only,
    sarima             = pred_sarima,
    seasonal_naive     = pred_snaive
  )
  
  common_rows <- complete.cases(prediction_df)
  
  if (sum(common_rows) < nrow(prediction_df)) {
    warning(
      split_label, ": only ",
      sum(common_rows), " of ", nrow(prediction_df),
      " validation weeks have predictions from all five models."
    )
  }
  
  prediction_eval <- prediction_df[common_rows, ]
  
  # -----------------------------
  # Collect results -- all models scored on the same N weeks
  # -----------------------------
  bind_rows(
    cbind(
      Split = split_label,
      Model = "Proposed quasi-Poisson (AR + Met Splines)",
      N = nrow(prediction_eval),
      calc_metrics(prediction_eval$observed, prediction_eval$proposed, train_sd)
    ),
    cbind(
      Split = split_label,
      Model = "Negative Binomial (AR + Met Splines)",
      N = nrow(prediction_eval),
      calc_metrics(prediction_eval$observed, prediction_eval$negative_binomial, train_sd)
    ),
    cbind(
      Split = split_label,
      Model = "AR(8) quasi-Poisson (No Meteorology)",
      N = nrow(prediction_eval),
      calc_metrics(prediction_eval$observed, prediction_eval$ar_only, train_sd)
    ),
    cbind(
      Split = split_label,
      Model = "Sequential SARIMA",
      N = nrow(prediction_eval),
      calc_metrics(prediction_eval$observed, prediction_eval$sarima, train_sd)
    ),
    cbind(
      Split = split_label,
      Model = "Seasonal Naive Baseline",
      N = nrow(prediction_eval),
      calc_metrics(prediction_eval$observed, prediction_eval$seasonal_naive, train_sd)
    )
  )
}

# ==============================================================================
# Rolling-origin CV benchmark comparison
# Same four folds used during lag selection
# ==============================================================================
cv_folds <- data.frame(
  fold        = 1:4,
  train_start = c(1, 1, 1, 1),
  train_end   = c(120, 132, 144, 156),
  valid_start = c(121, 133, 145, 157),
  valid_end   = c(132, 144, 156, 168)
)

cv_results_list <- vector("list", nrow(cv_folds))

for (k in seq_len(nrow(cv_folds))) {
  cat("\n=== Processing Fold", k, "===\n")
  cv_results_list[[k]] <- evaluate_one_split(
    train_start = cv_folds$train_start[k], train_end = cv_folds$train_end[k],
    test_start  = cv_folds$valid_start[k], test_end  = cv_folds$valid_end[k],
    split_label = paste0("Fold ", k)
  )
}

benchmark_cv_all <- bind_rows(cv_results_list)

benchmark_cv_summary <- benchmark_cv_all %>%
  group_by(Model) %>%
  summarise(
    mean_RMSE  = mean(RMSE, na.rm = TRUE),  sd_RMSE  = sd(RMSE, na.rm = TRUE),
    mean_MAE   = mean(MAE, na.rm = TRUE),   sd_MAE   = sd(MAE, na.rm = TRUE),
    mean_SRMSE = mean(SRMSE, na.rm = TRUE), sd_SRMSE = sd(SRMSE, na.rm = TRUE),
    mean_SMAE  = mean(SMAE, na.rm = TRUE),  sd_SMAE  = sd(SMAE, na.rm = TRUE),
    n_folds = sum(!is.na(RMSE)),
    .groups = "drop"
  ) %>%
  arrange(mean_SRMSE)

cat("\n--- Rolling-origin CV benchmark comparison: fold-level ---\n")
print(benchmark_cv_all)
cat("\n--- Rolling-origin CV benchmark comparison: mean +/- SD ---\n")
print(benchmark_cv_summary)

knitr::kable(benchmark_cv_summary, digits = 4,
             caption = "Rolling-origin cross-validation benchmark comparison")

write.csv(benchmark_cv_all,     file.path(out_dir, "benchmark_cv_fold_level_results.csv"), row.names = FALSE)
write.csv(benchmark_cv_summary, file.path(out_dir, "benchmark_cv_summary_results.csv"), row.names = FALSE)

cat("\nRolling-origin CV benchmark comparison completed. Outputs saved to:", out_dir, "\n")
