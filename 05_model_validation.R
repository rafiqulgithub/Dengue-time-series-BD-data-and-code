
# ============================================================

library(dplyr)
library(splines)
library(ggplot2)
library(tidyr)
library(broom)
library(knitr)
library(MASS)
library(forecast)
library(patchwork)

# -----------------------------
# Output directory -- change to your local folder
# -----------------------------
out_dir <- "~/Downloads/full_pipeline_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 0. Load data
# ============================================================
df <- read.csv("MyBDweeklyFinaldata.csv",stringsAsFactors = FALSE)

df <- df %>%
  dplyr::arrange(Week) %>%
  dplyr::mutate(
    time_index = row_number(),
    trend = time_index,
    sin1 = sin(2 * pi * time_index / 52),
    cos1 = cos(2 * pi * time_index / 52)
  )

final_ar_order   <- 8
final_temp_lag   <- 9
final_precip_lag <- 7
final_humid_lag  <- 8

# ============================================================
# Helper functions (shared across all sections)
# ============================================================
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

build_ns_basis <- function(train_x) {
  train_x <- train_x[is.finite(train_x)]
  knots <- as.numeric(quantile(train_x, probs = c(0.25, 0.50, 0.75), na.rm = TRUE))
  knots <- unique(knots)
  boundary_knots <- range(train_x, na.rm = TRUE)
  basis_train <- ns(train_x, knots = knots, Boundary.knots = boundary_knots)
  list(train_basis = basis_train, knots = knots, boundary_knots = boundary_knots,
       predict_fn = function(new_x) predict(basis_train, new_x))
}

add_spline_terms <- function(train_data, test_data = NULL) {
  temp_basis   <- build_ns_basis(train_data$Temp_lag)
  precip_basis <- build_ns_basis(train_data$Precip_lag)
  humid_basis  <- build_ns_basis(train_data$Humd_lag)
  
  train_temp <- as.data.frame(temp_basis$train_basis)
  train_prec <- as.data.frame(precip_basis$train_basis)
  train_humd <- as.data.frame(humid_basis$train_basis)
  
  colnames(train_temp) <- paste0("ns_temp_", seq_len(ncol(train_temp)))
  colnames(train_prec) <- paste0("ns_precip_", seq_len(ncol(train_prec)))
  colnames(train_humd) <- paste0("ns_humid_", seq_len(ncol(train_humd)))
  
  train_full <- bind_cols(train_data, train_temp, train_prec, train_humd)
  spline_terms <- c(colnames(train_temp), colnames(train_prec), colnames(train_humd))
  
  if (is.null(test_data)) {
    return(list(train_full = train_full, spline_terms = spline_terms,
                temp_basis = temp_basis, precip_basis = precip_basis, humid_basis = humid_basis))
  }
  
  test_temp <- as.data.frame(temp_basis$predict_fn(test_data$Temp_lag))
  test_prec <- as.data.frame(precip_basis$predict_fn(test_data$Precip_lag))
  test_humd <- as.data.frame(humid_basis$predict_fn(test_data$Humd_lag))
  
  colnames(test_temp) <- colnames(train_temp)
  colnames(test_prec) <- colnames(train_prec)
  colnames(test_humd) <- colnames(train_humd)
  
  test_full <- bind_cols(test_data, test_temp, test_prec, test_humd)
  
  list(train_full = train_full, test_full = test_full, spline_terms = spline_terms,
       temp_basis = temp_basis, precip_basis = precip_basis, humid_basis = humid_basis)
}

# ============================================================
# SECTION 1: Independent train-test validation (weeks 1-168 / 169-200)
# ============================================================
train_end_week  <- 168
test_start_week <- 169

max_lag <- max(final_ar_order, final_temp_lag, final_precip_lag, final_humid_lag)

train_raw <- df %>% dplyr::filter(time_index <= train_end_week)
test_raw  <- df %>% dplyr::filter(time_index >= train_end_week - max_lag + 1)

train_lagged <- make_model_data(train_raw, final_ar_order, final_temp_lag, final_precip_lag, final_humid_lag)
test_lagged  <- make_model_data(test_raw,  final_ar_order, final_temp_lag, final_precip_lag, final_humid_lag)

ar_terms <- paste0("lag_dengue_", seq_len(final_ar_order))
needed_cols <- c("Dengue_cases", ar_terms, "Temp_lag", "Precip_lag", "Humd_lag",
                 "trend", "sin1", "cos1")

train_model <- train_lagged %>% dplyr::filter(if_all(all_of(needed_cols), ~ !is.na(.)))
test_model <- test_lagged %>%
  dplyr::filter(if_all(all_of(needed_cols), ~ !is.na(.))) %>%
  dplyr::filter(time_index >= test_start_week)

cv_fold_max_index <- 168
overlap_n <- sum(test_model$time_index <= cv_fold_max_index)
cat("Overlap between test set and CV lag-selection folds:", overlap_n, "rows (should be 0)\n")
cat("Training weeks:", nrow(train_model), " | Independent test weeks:", nrow(test_model), "\n\n")

spline_obj   <- add_spline_terms(train_model, test_model)
train_full   <- spline_obj$train_full
test_full    <- spline_obj$test_full
spline_terms <- spline_obj$spline_terms

temp_terms   <- grep("^ns_temp_",   spline_terms, value = TRUE)
precip_terms <- grep("^ns_precip_", spline_terms, value = TRUE)
humid_terms  <- grep("^ns_humid_",  spline_terms, value = TRUE)

final_formula <- as.formula(
  paste("Dengue_cases ~", paste(ar_terms, collapse = " + "), "+",
        paste(spline_terms, collapse = " + "), "+", "trend + sin1 + cos1")
)

validation_fit <- glm(final_formula, family = quasipoisson(link = "log"), data = train_full)
cat("\n--- Training-period model summary ---\n")
print(summary(validation_fit))

# ------------------------------------------------------------
# 1a. Coefficient table WITH 95% confidence intervals
# CI uses the t-distribution with model residual df, consistent
# with the quasi-Poisson dispersion-adjusted standard errors
# (matches the approach used later in the partial-effect script).
# ------------------------------------------------------------
crit_val <- qt(0.975, df = validation_fit$df.residual)

coef_table <- broom::tidy(validation_fit) %>%
  dplyr::mutate(
    conf.low  = estimate - crit_val * std.error,
    conf.high = estimate + crit_val * std.error,
    estimate  = round(estimate, 4),
    std.error = round(std.error, 4),
    statistic = round(statistic, 3),
    conf.low  = round(conf.low, 4),
    conf.high = round(conf.high, 4),
    p.value   = ifelse(p.value < 0.001, "<0.001", as.character(round(p.value, 4)))
  ) %>%
  dplyr::rename(Term = term, Estimate = estimate, `Std. Error` = std.error,
                `t value` = statistic, `p-value` = p.value,
                `95% CI Lower` = conf.low, `95% CI Upper` = conf.high) %>%
  dplyr::select(Term, Estimate, `Std. Error`, `95% CI Lower`, `95% CI Upper`, `t value`, `p-value`)

cat("\n--- Coefficient table with 95% CI (training-period model) ---\n")
print(coef_table)
write.csv(coef_table, file.path(out_dir, "training_coefficient_table_with_CI.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 1b. Dispersion, deviance-based pseudo R-squared, squared
#     Pearson correlation (observed vs fitted, training)
# ------------------------------------------------------------
dispersion_parameter <- summary(validation_fit)$dispersion
pearson_dispersion <- sum(residuals(validation_fit, type = "pearson")^2) / validation_fit$df.residual
pseudo_R2_deviance <- 1 - validation_fit$deviance / validation_fit$null.deviance

train_pred <- predict(validation_fit, newdata = train_full, type = "response")
r2_pearson <- cor(train_full$Dengue_cases, train_pred)^2

model_summary_table <- data.frame(
  Metric = c("Observations used", "Residual degrees of freedom", "Null deviance",
             "Residual deviance", "Dispersion parameter", "Pearson dispersion",
             "Deviance-based pseudo R-squared", "Squared Pearson correlation (R2)"),
  Value = c(nobs(validation_fit), validation_fit$df.residual,
            round(validation_fit$null.deviance, 2), round(validation_fit$deviance, 2),
            round(dispersion_parameter, 3), round(pearson_dispersion, 3),
            round(pseudo_R2_deviance, 4), round(r2_pearson, 4))
)

cat("\n--- Model summary table (training period) ---\n")
print(model_summary_table)
write.csv(model_summary_table, file.path(out_dir, "training_fit_summary.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 1c. Joint significance (F-tests) for each predictor block
# ------------------------------------------------------------
formula_no_temp <- as.formula(
  paste("Dengue_cases ~", paste(ar_terms, collapse = " + "), "+",
        paste(precip_terms, collapse = " + "), "+",
        paste(humid_terms, collapse = " + "), "+", "trend + sin1 + cos1")
)
formula_no_precip <- as.formula(
  paste("Dengue_cases ~", paste(ar_terms, collapse = " + "), "+",
        paste(temp_terms, collapse = " + "), "+",
        paste(humid_terms, collapse = " + "), "+", "trend + sin1 + cos1")
)
formula_no_humid <- as.formula(
  paste("Dengue_cases ~", paste(ar_terms, collapse = " + "), "+",
        paste(temp_terms, collapse = " + "), "+",
        paste(precip_terms, collapse = " + "), "+", "trend + sin1 + cos1")
)
formula_no_ar <- as.formula(
  paste("Dengue_cases ~", paste(temp_terms, collapse = " + "), "+",
        paste(precip_terms, collapse = " + "), "+",
        paste(humid_terms, collapse = " + "), "+", "trend + sin1 + cos1")
)

model_no_temp   <- glm(formula_no_temp,   family = quasipoisson(link = "log"), data = train_full)
model_no_precip <- glm(formula_no_precip, family = quasipoisson(link = "log"), data = train_full)
model_no_humid  <- glm(formula_no_humid,  family = quasipoisson(link = "log"), data = train_full)
model_no_ar     <- glm(formula_no_ar,     family = quasipoisson(link = "log"), data = train_full)

cat("\n--- Joint significance: Temperature block ---\n")
print(anova(model_no_temp, validation_fit, test = "F"))
cat("\n--- Joint significance: Precipitation block ---\n")
print(anova(model_no_precip, validation_fit, test = "F"))
cat("\n--- Joint significance: Humidity block ---\n")
print(anova(model_no_humid, validation_fit, test = "F"))
cat("\n--- Joint significance: Autoregressive block ---\n")
print(anova(model_no_ar, validation_fit, test = "F"))

# ------------------------------------------------------------
# 1d. Out-of-sample validation performance
# ------------------------------------------------------------
test_pred <- predict(validation_fit, newdata = test_full, type = "response")
train_sd  <- sd(train_full$Dengue_cases, na.rm = TRUE)

validation_performance <- data.frame(
  Dataset = c("Training", "Independent testing"),
  N = c(nrow(train_full), nrow(test_full)),
  RMSE = c(rmse(train_full$Dengue_cases, train_pred), rmse(test_full$Dengue_cases, test_pred)),
  MAE  = c(mae(train_full$Dengue_cases, train_pred),  mae(test_full$Dengue_cases, test_pred))
) %>%
  dplyr::mutate(SRMSE = RMSE / train_sd, SMAE = MAE / train_sd) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4)))

cat("\n--- Out-of-sample validation performance ---\n")
print(validation_performance)
write.csv(validation_performance, file.path(out_dir, "validation_performance_table.csv"), row.names = FALSE)

obs_pred_df <- bind_rows(
  data.frame(time_index = train_full$time_index, observed = train_full$Dengue_cases,
             predicted = train_pred, period = "Training"),
  data.frame(time_index = test_full$time_index, observed = test_full$Dengue_cases,
             predicted = test_pred, period = "Independent testing")
)

obs_pred_plot <- ggplot(obs_pred_df, aes(x = time_index)) +
  geom_line(aes(y = observed, color = "Observed"), linewidth = 0.7) +
  geom_line(aes(y = predicted, color = "Predicted"), linewidth = 0.7, linetype = "dashed") +
  geom_vline(xintercept = train_end_week, linetype = "dotted") +
  labs(x = "Week", y = "Weekly dengue cases", color = "",
       title = "Observed versus predicted weekly dengue cases") +
  theme_minimal(base_size = 13)

print(obs_pred_plot)
ggsave(file.path(out_dir, "observed_vs_predicted_validation.png"), obs_pred_plot,
       width = 10, height = 5, dpi = 300)

# ============================================================
# SECTION 2: Partial-effect (exposure-response) plots
# ============================================================
make_effect_curve <- function(model, raw_x, basis_obj, spline_terms, variable_label,
                              reference_value = NULL, n_points = 200) {
  raw_x <- raw_x[is.finite(raw_x)]
  if (length(raw_x) == 0L) stop("No finite exposure values were supplied for ", variable_label, ".")
  n_points <- as.integer(n_points)
  
  missing_terms <- setdiff(spline_terms, names(coef(model)))
  if (length(missing_terms) > 0L) {
    stop("The following spline coefficients were not found in the fitted model: ",
         paste(missing_terms, collapse = ", "))
  }
  
  if (is.null(reference_value)) reference_value <- median(raw_x)
  if (!is.finite(reference_value)) stop("The reference value must be finite.")
  
  exposure_grid <- seq(min(raw_x), max(raw_x), length.out = n_points)
  basis_grid <- as.matrix(basis_obj$predict_fn(exposure_grid))
  basis_reference <- as.matrix(basis_obj$predict_fn(reference_value))
  
  if (ncol(basis_grid) != length(spline_terms)) {
    stop("The number of spline-basis columns does not match the number of supplied ",
         "spline coefficient names for ", variable_label, ".")
  }
  colnames(basis_grid) <- spline_terms
  
  beta <- coef(model)[spline_terms]
  covariance_matrix <- vcov(model)[spline_terms, spline_terms, drop = FALSE]
  
  centered_basis <- sweep(basis_grid, MARGIN = 2, STATS = as.numeric(basis_reference), FUN = "-")
  log_mean_ratio <- as.numeric(centered_basis %*% beta)
  effect_variance <- rowSums((centered_basis %*% covariance_matrix) * centered_basis)
  standard_error <- sqrt(pmax(effect_variance, 0))
  critical_value <- qt(0.975, df = model$df.residual)
  
  data.frame(
    variable = variable_label, exposure = exposure_grid, mean_ratio = exp(log_mean_ratio),
    lower_95 = exp(log_mean_ratio - critical_value * standard_error),
    upper_95 = exp(log_mean_ratio + critical_value * standard_error),
    reference_value = reference_value
  )
}

temperature_curve <- make_effect_curve(
  validation_fit, train_full$Temp_lag, spline_obj$temp_basis, temp_terms,
  paste0("Temperature (lag ", final_temp_lag, " weeks)")
)
precipitation_curve <- make_effect_curve(
  validation_fit, train_full$Precip_lag, spline_obj$precip_basis, precip_terms,
  paste0("Precipitation (lag ", final_precip_lag, " weeks)")
)
humidity_curve <- make_effect_curve(
  validation_fit, train_full$Humd_lag, spline_obj$humid_basis, humid_terms,
  paste0("Relative humidity (lag ", final_humid_lag, " weeks)")
)

plot_single_effect <- function(curve_data, rug_values, x_label, plot_title) {
  reference_value <- unique(curve_data$reference_value)
  ggplot(curve_data, aes(x = exposure, y = mean_ratio)) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.25) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6) +
    geom_vline(xintercept = reference_value, linetype = "dotted", linewidth = 0.6) +
    geom_rug(data = data.frame(exposure = rug_values), aes(x = exposure),
             inherit.aes = FALSE, sides = "b", alpha = 0.30) +
    labs(x = x_label, y = "Fitted mean ratio of weekly dengue incidence", title = plot_title) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
}

temperature_plot <- plot_single_effect(
  temperature_curve, train_full$Temp_lag,
  paste0("Weekly mean temperature at lag ", final_temp_lag, " (\u00B0C)"), "(a) Temperature"
)
precipitation_plot <- plot_single_effect(
  precipitation_curve, train_full$Precip_lag,
  paste0("Weekly cumulative precipitation at lag ", final_precip_lag, " (mm)"), "(b) Precipitation"
)
humidity_plot <- plot_single_effect(
  humidity_curve, train_full$Humd_lag,
  paste0("Weekly mean relative humidity at lag ", final_humid_lag, " (%)"), "(c) Relative humidity"
)

combined_partial_effect_plot <-
  temperature_plot + (precipitation_plot + ylab(NULL)) + (humidity_plot + ylab(NULL)) +
  plot_layout(nrow = 1)

print(combined_partial_effect_plot)
ggsave(file.path(out_dir, "partial_effect_plots_combined.png"), combined_partial_effect_plot,
       width = 15, height = 5, dpi = 600, bg = "white")

effect_curves_all <- bind_rows(temperature_curve, precipitation_curve, humidity_curve)
write.csv(effect_curves_all, file.path(out_dir, "partial_effect_curve_estimates.csv"), row.names = FALSE)

# ============================================================
# SECTION 3: Residual diagnostics (training period)
# ============================================================
diag_df <- data.frame(
  week     = train_full$time_index,
  pearson  = residuals(validation_fit, type = "pearson"),
  deviance = residuals(validation_fit, type = "deviance"),
  observed = train_full$Dengue_cases,
  fitted   = fitted(validation_fit)
)

test_lags <- c(6, 12, 20)
ljung_results <- lapply(test_lags, function(L) {
  test <- Box.test(diag_df$pearson, lag = L, type = "Ljung-Box")
  data.frame(Lag = L, Statistic = unname(test$statistic),
             df = unname(test$parameter), p_value = test$p.value)
})
ljung_table <- do.call(rbind, ljung_results)
cat("\n--- Ljung-Box test: training-period residuals ---\n")
print(ljung_table)
write.csv(ljung_table, file.path(out_dir, "ljung_box_residual_test.csv"), row.names = FALSE)

plot_a <- ggplot(diag_df, aes(x = week, y = pearson)) +
  geom_line(color = "black", linewidth = 0.75) +
  geom_hline(yintercept = 0, color = "red3", linewidth = 0.9) +
  labs(title = "Residual Sequence", x = "Week", y = "Pearson Residuals") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

plot_b <- ggplot(diag_df, aes(x = deviance)) +
  geom_histogram(aes(y = after_stat(density)), bins = 18,
                 fill = "grey85", color = "grey35", linewidth = 0.8) +
  geom_density(color = "#005AB5", linewidth = 1.4) +
  labs(title = "Residual Distribution", x = "Deviance Residuals", y = "Density") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

plot_c <- ggplot(diag_df, aes(sample = deviance)) +
  stat_qq(size = 2.2) +
  stat_qq_line(color = "red3", linewidth = 1) +
  labs(title = "Q-Q Plot", x = "Theoretical Quantiles", y = "Deviance Residuals") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

plot_d <- ggplot(diag_df, aes(x = fitted, y = pearson)) +
  geom_point(size = 2.3, alpha = 0.70, shape = 16, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.9) +
  geom_smooth(method = "loess", se = FALSE, color = "red3", linewidth = 1.2) +
  labs(title = "Residuals vs Fitted", x = "Fitted Values", y = "Pearson Residuals") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

acf_vals <- acf(diag_df$pearson, plot = FALSE, lag.max = 20)
acf_df <- data.frame(lag = as.numeric(acf_vals$lag), acf = as.numeric(acf_vals$acf))
ci_line <- qnorm(0.975) / sqrt(nrow(diag_df))

plot_e <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
  geom_hline(yintercept = c(-ci_line, ci_line), linetype = "dashed",
             linewidth = 0.9, color = "#005AB5") +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 1.0, color = "black") +
  geom_point(size = 2.2, color = "black") +
  labs(title = "Residual ACF", x = "Lag (weeks)", y = "Autocorrelation") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

diagnostic_panel_5 <- (plot_a | plot_b | plot_c) /
  (plot_d | plot_e | plot_spacer()) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 13))

print(diagnostic_panel_5)
ggsave(file.path(out_dir, "residual_diagnostics_panel_5.png"), diagnostic_panel_5,
       width = 13, height = 8, dpi = 300)


# ============================================================
# SECTION 4: Refit final model on complete data + 7-week forecast
# ============================================================
full_lagged <- make_model_data(df, final_ar_order, final_temp_lag, final_precip_lag, final_humid_lag)
full_model  <- full_lagged %>% dplyr::filter(if_all(all_of(needed_cols), ~ !is.na(.)))

full_spline_obj   <- add_spline_terms(full_model)
full_data         <- full_spline_obj$train_full
full_spline_terms <- full_spline_obj$spline_terms

full_formula <- as.formula(
  paste("Dengue_cases ~", paste(ar_terms, collapse = " + "), "+",
        paste(full_spline_terms, collapse = " + "), "+", "trend + sin1 + cos1")
)

final_full_fit <- glm(full_formula, family = quasipoisson(link = "log"), data = full_data)
cat("\n--- Final model (refit on complete data) summary ---\n")
print(summary(final_full_fit))

forecast_horizon <- min(final_temp_lag, final_precip_lag, final_humid_lag)
last_observed_time <- max(df$time_index)

predictiondata <- data.frame(
  forecast_week = seq_len(forecast_horizon),
  time_index    = last_observed_time + seq_len(forecast_horizon),
  Temp_lag      = df$Temperature[last_observed_time + seq_len(forecast_horizon) - final_temp_lag],
  Precip_lag    = df$Precipitation[last_observed_time + seq_len(forecast_horizon) - final_precip_lag],
  Humd_lag      = df$Humidity[last_observed_time + seq_len(forecast_horizon) - final_humid_lag]
)
predictiondata$trend <- predictiondata$time_index
predictiondata$sin1  <- sin(2 * pi * predictiondata$time_index / 52)
predictiondata$cos1  <- cos(2 * pi * predictiondata$time_index / 52)

cat("\n--- Lag-aligned meteorological predictors for forecasting ---\n")
print(predictiondata)

dengue_history <- tail(df$Dengue_cases, final_ar_order)
forecast_values <- numeric(forecast_horizon)

full_temp_terms   <- grep("^ns_temp_",   full_spline_terms, value = TRUE)
full_precip_terms <- grep("^ns_precip_", full_spline_terms, value = TRUE)
full_humid_terms  <- grep("^ns_humid_",  full_spline_terms, value = TRUE)

for (h in seq_len(forecast_horizon)) {
  newrow <- data.frame(
    time_index = predictiondata$time_index[h], trend = predictiondata$trend[h],
    sin1 = predictiondata$sin1[h], cos1 = predictiondata$cos1[h],
    Temp_lag = predictiondata$Temp_lag[h], Precip_lag = predictiondata$Precip_lag[h],
    Humd_lag = predictiondata$Humd_lag[h]
  )
  for (i in seq_len(final_ar_order)) {
    newrow[[paste0("lag_dengue_", i)]] <- log(dengue_history[final_ar_order + 1 - i] + 1)
  }
  
  temp_spline   <- as.data.frame(full_spline_obj$temp_basis$predict_fn(newrow$Temp_lag))
  precip_spline <- as.data.frame(full_spline_obj$precip_basis$predict_fn(newrow$Precip_lag))
  humid_spline  <- as.data.frame(full_spline_obj$humid_basis$predict_fn(newrow$Humd_lag))
  
  colnames(temp_spline)   <- full_temp_terms
  colnames(precip_spline) <- full_precip_terms
  colnames(humid_spline)  <- full_humid_terms
  
  newrow_full <- dplyr::bind_cols(newrow, temp_spline, precip_spline, humid_spline)
  
  forecast_values[h] <- as.numeric(predict(final_full_fit, newdata = newrow_full, type = "response"))
  dengue_history <- c(dengue_history[-1], forecast_values[h])
}

forecast_results <- data.frame(
  forecast_week   = seq_len(forecast_horizon),
  time_index      = predictiondata$time_index,
  predicted_cases = round(forecast_values, 0)
)

observed_seven_weeks <- c(6177, 6552, 5827, 4301, 3638, 2866, 1656)
if (length(observed_seven_weeks) != forecast_horizon) {
  stop("observed_seven_weeks length (", length(observed_seven_weeks),
       ") does not match forecast_horizon (", forecast_horizon, ")")
}

forecast_results$observed_cases <- observed_seven_weeks
forecast_results$Year_week <- c("2025 W44", "2025 W45", "2025 W46", "2025 W47",
                                "2025 W48", "2025 W49", "2025 W50")

cat("\n--- Forecast versus observed cases ---\n")
print(forecast_results)
write.csv(forecast_results, file.path(out_dir, "seven_week_forecast_results.csv"), row.names = FALSE)

forecast_rmse <- rmse(forecast_results$observed_cases, forecast_results$predicted_cases)
forecast_mae  <- mae(forecast_results$observed_cases,  forecast_results$predicted_cases)
forecast_mape <- mean(abs((forecast_results$observed_cases - forecast_results$predicted_cases)
                          / forecast_results$observed_cases), na.rm = TRUE) * 100

forecast_performance <- data.frame(RMSE = forecast_rmse, MAE = forecast_mae, MAPE = forecast_mape) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("\n--- Seven-week forecast performance ---\n")
print(forecast_performance)
write.csv(forecast_performance, file.path(out_dir, "seven_week_forecast_performance.csv"), row.names = FALSE)

forecast_plot_data <- forecast_results %>%
  dplyr::select(Year_week, predicted_cases, observed_cases) %>%
  dplyr::rename(Forecast = predicted_cases, Observed = observed_cases) %>%
  tidyr::pivot_longer(cols = c(Forecast, Observed), names_to = "Type", values_to = "Cases")

forecast_plot_data$Year_week <- factor(forecast_plot_data$Year_week, levels = unique(forecast_results$Year_week))
forecast_plot_data$Type <- factor(forecast_plot_data$Type, levels = c("Forecast", "Observed"))

seven_week_forecast_plot <- ggplot(forecast_plot_data,
                                   aes(x = Year_week, y = Cases, color = Type, linetype = Type, group = Type)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Forecast" = "blue", "Observed" = "red")) +
  scale_linetype_manual(values = c("Forecast" = "solid", "Observed" = "dashed")) +
  labs(x = "Year-week", y = "Weekly dengue cases", color = NULL, linetype = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 1)
  )

print(seven_week_forecast_plot)
ggsave(file.path(out_dir, "seven_week_ahead_forecast.png"), seven_week_forecast_plot,
       width = 10, height = 6, dpi = 300)

cat("\n============================================================\n")
cat("Pipeline completed. All outputs saved to:", out_dir, "\n")
cat("============================================================\n")