
library(dplyr)
library(splines)

df <-   read.csv("MyBDweeklyFinaldata.csv",stringsAsFactors = FALSE)
head(df)

df <- df %>%
  dplyr::arrange(Week) %>%
  dplyr::mutate(
    time_index = dplyr::row_number(),
    trend = time_index,
    sin1 = sin(2 * pi * time_index / 52),
    cos1 = cos(2 * pi * time_index / 52)
  )

# -----------------------------
# Helper functions
# -----------------------------
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae  <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)

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

# Build a natural spline basis from TRAINING data knots, then apply
# the identical basis (same knots) to validation data. This avoids
# leaking validation-period quantiles into the spline construction.
build_ns_basis <- function(train_x, valid_x) {
  train_x <- train_x[is.finite(train_x)]
  kn <- stats::quantile(train_x, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
  bk <- range(train_x, na.rm = TRUE)
  # ns() needs distinct knots; guard against degenerate cases
  kn <- unique(kn)
  basis_train <- ns(train_x, knots = kn, Boundary.knots = bk)
  list(
    train = basis_train,
    predict_fn = function(newx) {
      predict(basis_train, newx)
    }
  )
}

# -----------------------------
# Candidate model grid
# (kept as your DLNM-informed ranges)
# -----------------------------
candidate_grid <- expand.grid(
  ar_order   = 2:12,
  temp_lag   = 8:17,
  precip_lag = 5:12,
  humid_lag  = 2:8
)

# -----------------------------
# Forward-chaining 4 folds
# -----------------------------
cv_folds <- data.frame(
  fold        = 1:4,
  train_start = c(1, 1, 1, 1),
  train_end   = c(120, 132, 144, 156),
  valid_start = c(121, 133, 145, 157),
  valid_end   = c(132, 144, 156, 168)

)

# -----------------------------
# Main forward-chaining CV loop
# -----------------------------
cv_results <- vector("list", nrow(candidate_grid))

for (g in seq_len(nrow(candidate_grid))) {
  
  ar_order   <- candidate_grid$ar_order[g]
  temp_lag   <- candidate_grid$temp_lag[g]
  precip_lag <- candidate_grid$precip_lag[g]
  humid_lag  <- candidate_grid$humid_lag[g]
  
  max_lag <- max(ar_order, temp_lag, precip_lag, humid_lag)
  
  fold_results <- vector("list", nrow(cv_folds))
  
  for (k in seq_len(nrow(cv_folds))) {
    
    train_start <- cv_folds$train_start[k]
    train_end   <- cv_folds$train_end[k]
    valid_start <- cv_folds$valid_start[k]
    valid_end   <- cv_folds$valid_end[k]
    
    train_raw <- df[train_start:train_end, ]
    valid_raw <- df[(train_end - max_lag + 1):valid_end, ]
    
    train_lagged <- make_model_data(train_raw, ar_order, temp_lag, precip_lag, humid_lag)
    valid_lagged <- make_model_data(valid_raw, ar_order, temp_lag, precip_lag, humid_lag)
    
    ar_terms <- paste0("lag_dengue_", seq_len(ar_order))
    
    needed_cols <- c("Dengue_cases", ar_terms, "Temp_lag", "Precip_lag", "Humd_lag",
                     "trend", "sin1", "cos1")
    
    train_model <- train_lagged %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(needed_cols), ~ !is.na(.)))
    
    valid_model <- valid_lagged %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(needed_cols), ~ !is.na(.))) %>%
      dplyr::filter(time_index >= valid_start & time_index <= valid_end)
    
    if (nrow(train_model) < 20 || nrow(valid_model) < 1) {
      fold_results[[k]] <- NULL
      next
    }
    
    # --- build spline bases from TRAINING data only, apply to both ---
    temp_basis   <- build_ns_basis(train_model$Temp_lag,   valid_model$Temp_lag)
    precip_basis <- build_ns_basis(train_model$Precip_lag, valid_model$Precip_lag)
    humid_basis  <- build_ns_basis(train_model$Humd_lag,   valid_model$Humd_lag)
    
    train_spline_df <- as.data.frame(temp_basis$train)
    colnames(train_spline_df) <- paste0("ns_temp_", seq_len(ncol(train_spline_df)))
    train_spline_df2 <- as.data.frame(predict(precip_basis$train, train_model$Precip_lag))
    colnames(train_spline_df2) <- paste0("ns_precip_", seq_len(ncol(train_spline_df2)))
    train_spline_df3 <- as.data.frame(predict(humid_basis$train, train_model$Humd_lag))
    colnames(train_spline_df3) <- paste0("ns_humid_", seq_len(ncol(train_spline_df3)))
    
    valid_spline_df  <- as.data.frame(temp_basis$predict_fn(valid_model$Temp_lag))
    colnames(valid_spline_df) <- paste0("ns_temp_", seq_len(ncol(valid_spline_df)))
    valid_spline_df2 <- as.data.frame(precip_basis$predict_fn(valid_model$Precip_lag))
    colnames(valid_spline_df2) <- paste0("ns_precip_", seq_len(ncol(valid_spline_df2)))
    valid_spline_df3 <- as.data.frame(humid_basis$predict_fn(valid_model$Humd_lag))
    colnames(valid_spline_df3) <- paste0("ns_humid_", seq_len(ncol(valid_spline_df3)))
    
    train_full <- dplyr::bind_cols(train_model, train_spline_df, train_spline_df2, train_spline_df3)
    valid_full <- dplyr::bind_cols(valid_model, valid_spline_df, valid_spline_df2, valid_spline_df3)
    
    spline_terms <- c(colnames(train_spline_df), colnames(train_spline_df2), colnames(train_spline_df3))
    
    model_formula <- as.formula(
      paste("Dengue_cases ~",
            paste(ar_terms, collapse = " + "), "+",
            paste(spline_terms, collapse = " + "), "+",
            "trend + sin1 + cos1")
    )
    
    fit <- tryCatch(
      glm(model_formula, family = quasipoisson(link = "log"), data = train_full),
      error = function(e) NULL
    )
    
    if (is.null(fit)) { fold_results[[k]] <- NULL; next }
    
    valid_full$predicted <- suppressWarnings(
      predict(fit, newdata = valid_full, type = "response")
    )
    
    train_sd <- sd(train_model$Dengue_cases, na.rm = TRUE)
    
    
    fold_results[[k]] <- data.frame(
      ar_order = ar_order, temp_lag = temp_lag, precip_lag = precip_lag, humid_lag = humid_lag,
      fold = k, train_end = train_end, valid_start = valid_start, valid_end = valid_end,
      n_train = nrow(train_model),
      RMSE = rmse(valid_full$Dengue_cases, valid_full$predicted),
      MAE  = mae(valid_full$Dengue_cases, valid_full$predicted),
      standardized_RMSE = rmse(valid_full$Dengue_cases, valid_full$predicted) /  train_sd,
      standardized_MAE  = mae(valid_full$Dengue_cases, valid_full$predicted) /  train_sd
    )
  }
  
  cv_results[[g]] <- dplyr::bind_rows(fold_results)
}

cv_all_results <- dplyr::bind_rows(cv_results)

# -----------------------------
# Summarize performance across folds
# -----------------------------
cv_summary <- cv_all_results %>%
  dplyr::group_by(ar_order, temp_lag, precip_lag, humid_lag) %>%
  dplyr::summarise(
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    mean_MAE  = mean(MAE, na.rm = TRUE),
    mean_standardized_RMSE = mean(standardized_RMSE, na.rm = TRUE),
    mean_standardized_MAE  = mean(standardized_MAE, na.rm = TRUE),
    sd_standardized_RMSE   = sd(standardized_RMSE, na.rm = TRUE),
    min_n_train = min(n_train, na.rm = TRUE),
    n_folds = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_folds == nrow(cv_folds)) %>%   # require all folds succeeded
  dplyr::mutate(lag_sum = ar_order + temp_lag + precip_lag + humid_lag) %>%
  dplyr::arrange(mean_standardized_RMSE)

# -----------------------------
# Tolerance / parsimony rule ("1-SE style")
# Among candidates within `tol` of the best mean standardized RMSE,
# pick the most parsimonious (smallest ar_order, then smallest lag_sum)
# -----------------------------
tol <- 0.02   # 2% tolerance band -- adjust as appropriate and justify in methods

best_value <- min(cv_summary$mean_standardized_RMSE)
within_tol <- cv_summary %>%
  dplyr::filter(mean_standardized_RMSE <= best_value * (1 + tol)) %>%
  dplyr::arrange(ar_order, lag_sum, mean_standardized_RMSE)

selected_model <- within_tol %>% dplyr::slice(1)

cat("\n--- Raw best (minimum mean standardized RMSE) ---\n")
print(cv_summary %>% dplyr::slice(1))

cat("\n--- Candidates within", tol * 100, "% tolerance of best (parsimony pool) ---\n")
print(within_tol)

cat("\n--- Selected model (parsimony rule applied) ---\n")
print(selected_model)

# -----------------------------
# Stability diagnostic: does the selected lag combo perform
# consistently across all 4 folds, or is it driven by one fold?
# -----------------------------
selected_fold_detail <- cv_all_results %>%
  dplyr::filter(
    ar_order   == selected_model$ar_order,
    temp_lag   == selected_model$temp_lag,
    precip_lag == selected_model$precip_lag,
    humid_lag  == selected_model$humid_lag
  )

cat("\n--- Per-fold performance of selected model (check for instability) ---\n")
print(selected_fold_detail)

# Save full results for inspection / sensitivity checks
write.csv(cv_all_results, "/home/claude/cv_all_results.csv", row.names = FALSE)
write.csv(cv_summary,     "/home/claude/cv_summary.csv", row.names = FALSE)