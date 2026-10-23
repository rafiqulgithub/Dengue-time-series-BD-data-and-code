
rm(list = ls())
library(dplyr)
library(splines)
library(ggplot2)
library(tidyr)
library(broom)
library(knitr)
library(patchwork)

# -----------------------------
# Output directory -- change to your local folder
# -----------------------------
out_dir <- "~/Downloads/full_model_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 0. Load data
# ============================================================
df <- read.csv("MyBDweeklyFinaldata.csv", stringsAsFactors = FALSE)

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
# SECTION 4: RECURSIVE FORECASTING + EMPIRICAL PI CALIBRATION
# 26 evaluation origins (weeks 168-193), each with 7 recursive steps
# ============================================================
# Signed errors = observed - predicted.
# Calibration uses development origins 100-161 only.
# At every origin, refit the model using data available through that origin.
# For horizons > 1, use previous forecast means in the AR history.
# Weather inputs are lagged so that they are observed by the forecast origin.

forecast_horizon <- 7L
calibration_origins <- 100:161
evaluation_origins <- 168:193
origin_grid <- evaluation_origins

stopifnot(
  nrow(df) == 200L,
  identical(as.integer(df$time_index), seq_len(nrow(df))),
  max(calibration_origins) + forecast_horizon <= train_end_week,
  max(evaluation_origins) + forecast_horizon <= nrow(df),
  all(is.finite(df$Dengue_cases)),
  all(df$Dengue_cases >= 0),
  all(is.finite(as.matrix(df[, c("Temperature", "Precipitation", "Humidity")])))
)

# Fit the full quasi-Poisson model using observations up to one origin.
fit_origin <- function(origin) {
  historical <- df[seq_len(origin), , drop = FALSE]
  lagged <- make_model_data(
    historical, final_ar_order,
    final_temp_lag, final_precip_lag, final_humid_lag
  )
  eligible <- lagged[complete.cases(lagged[, needed_cols, drop = FALSE]), , drop = FALSE]
  if (nrow(eligible) < 50L) {
    stop("Too few complete training rows at origin ", origin)
  }

  bases <- add_spline_terms(eligible)
  form <- stats::reformulate(
    c(ar_terms, bases$spline_terms, "trend", "sin1", "cos1"),
    response = "Dengue_cases"
  )
  fit <- stats::glm(
    form, family = stats::quasipoisson(link = "log"),
    data = bases$train_full
  )
  if (anyNA(stats::coef(fit)) || any(!is.finite(stats::vcov(fit)))) {
    stop("Non-estimable coefficients at origin ", origin)
  }
  if (!isTRUE(fit$converged)) stop("GLM did not converge at origin ", origin)
  list(fit = fit, bases = bases)
}

# Construct a one-row prediction data frame from historical weather and
# the current AR history, ordered oldest -> newest.
future_design <- function(obj, origin, h, history8) {
  target <- origin + h
  weather_index <- c(
    target - final_temp_lag,
    target - final_precip_lag,
    target - final_humid_lag
  )
  if (any(weather_index < 1L) || any(weather_index > origin)) {
    stop("Unavailable weather / possible leakage at origin ", origin)
  }
  if (length(history8) != final_ar_order) {
    stop("AR history length is not equal to final_ar_order")
  }

  newdata <- data.frame(
    trend = target,
    sin1 = sin(2 * pi * target / 52),
    cos1 = cos(2 * pi * target / 52),
    Temp_lag = df$Temperature[weather_index[1]],
    Precip_lag = df$Precipitation[weather_index[2]],
    Humd_lag = df$Humidity[weather_index[3]]
  )
  for (i in seq_len(final_ar_order)) {
    newdata[[paste0("lag_dengue_", i)]] <-
      log1p(history8[final_ar_order + 1L - i])
  }

  spl <- obj$bases
  temp_terms <- grep("^ns_temp_", spl$spline_terms, value = TRUE)
  precip_terms <- grep("^ns_precip_", spl$spline_terms, value = TRUE)
  humid_terms <- grep("^ns_humid_", spl$spline_terms, value = TRUE)

  b_temp <- as.data.frame(spl$temp_basis$predict_fn(newdata$Temp_lag))
  b_precip <- as.data.frame(spl$precip_basis$predict_fn(newdata$Precip_lag))
  b_humid <- as.data.frame(spl$humid_basis$predict_fn(newdata$Humd_lag))
  names(b_temp) <- temp_terms
  names(b_precip) <- precip_terms
  names(b_humid) <- humid_terms
  dplyr::bind_cols(newdata, b_temp, b_precip, b_humid)
}

# Recursive point forecasts. Observations are attached for scoring only;
# they are never fed into the forecast recursion.
recursive_origin <- function(origin, observed = NULL) {
  obj <- fit_origin(origin)
  history8 <- tail(df$Dengue_cases[seq_len(origin)], final_ar_order)
  point <- numeric(forecast_horizon)

  for (h in seq_len(forecast_horizon)) {
    newdata <- future_design(obj, origin, h, history8)
    mu <- as.numeric(stats::predict(obj$fit, newdata = newdata, type = "response"))
    if (length(mu) != 1L || !is.finite(mu) || mu < 0) {
      stop("Invalid forecast at origin ", origin, ", horizon ", h)
    }
    point[h] <- mu
    history8 <- c(history8[-1], mu)
  }

  targets <- origin + seq_len(forecast_horizon)
  if (is.null(observed)) {
    observed <- if (max(targets) <= nrow(df)) {
      df$Dengue_cases[targets]
    } else {
      rep(NA_real_, forecast_horizon)
    }
  }
  if (length(observed) != forecast_horizon) {
    stop("Observed vector must have length forecast_horizon")
  }

  data.frame(
    model = "full",
    origin = origin,
    horizon = seq_len(forecast_horizon),
    target_week = targets,
    observed = as.numeric(observed),
    predicted = point
  )
}

# ------------------------------------------------------------
# 4A. Historical forecast-error calibration: origins 100-161
# ------------------------------------------------------------
cat("\n--- Calibrating empirical forecast errors (development only) ---\n")
calibration_results <- dplyr::bind_rows(lapply(calibration_origins, function(origin) {
  cat("Calibration origin", origin, "of", max(calibration_origins), "\n")
  out <- recursive_origin(origin)
  out$signed_error <- out$observed - out$predicted
  out
}))
stopifnot(all(calibration_results$target_week <= train_end_week))

# Horizon-specific empirical error quantiles. type = 1 matches script 06.
# These are empirical forecast-error intervals, not formal guaranteed 95% PIs.
error_quantiles <- calibration_results %>%
  dplyr::group_by(horizon) %>%
  dplyr::summarise(
    n_calibration = dplyr::n(),
    q025 = unname(stats::quantile(signed_error, 0.025, na.rm = TRUE, type = 1)),
    q975 = unname(stats::quantile(signed_error, 0.975, na.rm = TRUE, type = 1)),
    calibration_MAE = mean(abs(signed_error), na.rm = TRUE),
    .groups = "drop"
  )
print(error_quantiles)
write.csv(calibration_results,
          file.path(out_dir, "development_calibration_forecast_errors.csv"),
          row.names = FALSE)
write.csv(error_quantiles,
          file.path(out_dir, "empirical_horizon_error_quantiles.csv"),
          row.names = FALSE)

# Apply the calibration quantiles unchanged to evaluation origins and final forecast.
add_empirical_PI <- function(forecasts) {
  out <- dplyr::left_join(forecasts, error_quantiles, by = "horizon")
  if (anyNA(out$q025) || anyNA(out$q975)) stop("Missing horizon quantiles")
  out %>% dplyr::mutate(
    lower95 = pmax(0, predicted + q025),
    upper95 = pmax(0, predicted + q975)
  )
}

# ------------------------------------------------------------
# 4B. Primary rolling-origin evaluation: 26 origins x 7 horizons
# ------------------------------------------------------------
cat("\n--- 26-origin recursive evaluation ---\n")
backtest_results <- dplyr::bind_rows(lapply(evaluation_origins, function(origin) {
  cat("Evaluation origin", origin, "of", max(evaluation_origins), "\n")
  recursive_origin(origin)
})) %>%
  add_empirical_PI() %>%
  dplyr::mutate(
    error = observed - predicted,
    absolute_error = abs(error),
    squared_error = error^2,
    ape = dplyr::if_else(observed > 0, 100 * absolute_error / observed, NA_real_),
    covered95 = observed >= lower95 & observed <= upper95,
    interval_width = upper95 - lower95
  )

stopifnot(
  nrow(backtest_results) == length(evaluation_origins) * forecast_horizon,
  dplyr::n_distinct(backtest_results$origin) == 26L,
  all(!is.na(backtest_results$observed)),
  all(backtest_results$lower95 <= backtest_results$upper95)
)

summarise_forecast_performance <- function(dat) {
  dat %>% dplyr::summarise(
    n = dplyr::n(),
    RMSE = sqrt(mean(squared_error, na.rm = TRUE)),
    MAE = mean(absolute_error, na.rm = TRUE),
    MAPE = mean(ape, na.rm = TRUE),
    Mean_Error = mean(error, na.rm = TRUE),
    Covered = sum(covered95, na.rm = TRUE),
    Coverage_percent = 100 * mean(covered95, na.rm = TRUE),
    Mean_PI_width = mean(interval_width, na.rm = TRUE),
    .groups = "drop"
  )
}

performance_by_origin <- backtest_results %>%
  dplyr::group_by(origin) %>% summarise_forecast_performance()
performance_by_horizon <- backtest_results %>%
  dplyr::group_by(horizon) %>% summarise_forecast_performance()
performance_overall <- backtest_results %>%
  summarise_forecast_performance()

# Keep the legacy horizon_summary object used by the existing plotting code below.
horizon_summary <- performance_by_horizon %>%
  dplyr::transmute(
    horizon,
    N_origins = n,
    RMSE,
    MAE,
    MAPE,
    Coverage95 = Coverage_percent / 100,
    Mean_interval_width = Mean_PI_width
  )

#cat("\n--- 26-origin performance by origin ---\n")
#print(performance_by_origin, n = Inf)
#cat("\n--- 26-origin performance by horizon ---\n")
#print(performance_by_horizon, n = Inf)
#cat("\n--- 26-origin overall performance ---\n")
#print(performance_overall)
#print(horizon_summary, n = Inf)

write.csv(backtest_results,
          file.path(out_dir, "rolling_origin_recursive_predictions.csv"),
          row.names = FALSE)
write.csv(performance_by_origin,
          file.path(out_dir, "rolling_origin_performance_by_origin.csv"),
          row.names = FALSE)
write.csv(performance_by_horizon,
          file.path(out_dir, "rolling_origin_performance_by_horizon.csv"),
          row.names = FALSE)
write.csv(performance_overall,
          file.path(out_dir, "rolling_origin_performance_overall.csv"),
          row.names = FALSE)
write.csv(horizon_summary,
          file.path(out_dir, "rolling_origin_horizon_summary.csv"),
          row.names = FALSE)

p_rmse <- ggplot(horizon_summary, aes(horizon, RMSE)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  scale_x_continuous(breaks = seq_len(forecast_horizon)) +
  labs(x = "Forecast horizon (weeks)", y = "RMSE (cases)",
       title = "Recursive forecasting error") + theme_bw(base_size = 12)

p_coverage <- ggplot(horizon_summary, aes(horizon, Coverage95)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_x_continuous(breaks = seq_len(forecast_horizon)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Forecast horizon (weeks)", y = "Empirical PI coverage",
       title = "Prediction interval coverage") + theme_bw(base_size = 12)

#print(p_rmse + p_coverage)
ggsave(file.path(out_dir, "rolling_origin_accuracy_coverage.png"),
       p_rmse + p_coverage, width = 11, height = 4.5, dpi = 300)

# ============================================================
# SECTION 5: FINAL FORECAST AFTER REFITTING ON ALL 200 WEEKS
# ============================================================
# The observed future counts are used for evaluation ONLY.
# Interval quantiles remain those from calibration origins 100-161.

observed_seven_weeks <- c(6177, 6552, 5827, 4301, 3638, 2866, 1656)
stopifnot(length(observed_seven_weeks) == forecast_horizon)
# Verify these correspond to Nov 1-7, Nov 8-14, ..., Dec 13-19, 2025.

final_forecast <- recursive_origin(nrow(df), observed = observed_seven_weeks)
final_forecast <- add_empirical_PI(final_forecast)
final_forecast$study_week <- paste0("Week ", final_forecast$target_week)
final_forecast$start_date <- as.Date("2025-11-01") +
  7L * (final_forecast$horizon - 1L)
final_forecast$end_date <- final_forecast$start_date + 6L
final_forecast$period <- paste(format(final_forecast$start_date, "%b %d"),
                               format(final_forecast$end_date, "%b %d"), sep = "–")

final_performance <- data.frame(
  RMSE = rmse(final_forecast$observed, final_forecast$predicted),
  MAE = mae(final_forecast$observed, final_forecast$predicted),
  MAPE = 100 * mean(abs(final_forecast$observed - final_forecast$predicted) /
                      final_forecast$observed),
  Coverage95 = mean(final_forecast$observed >= final_forecast$lower95 &
                      final_forecast$observed <= final_forecast$upper95),
  Mean_interval_width = mean(final_forecast$upper95 - final_forecast$lower95)
)
cat("\n--- Final forecast with empirical prediction intervals ---\n")
print(final_forecast)
print(final_performance)
write.csv(final_forecast, file.path(out_dir, "seven_week_forecast_with_PI.csv"),
          row.names = FALSE)
write.csv(final_performance,
          file.path(out_dir, "seven_week_forecast_performance.csv"),
          row.names = FALSE)

final_plot <- ggplot(final_forecast, aes(x = horizon)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95),
              fill = "steelblue", alpha = 0.20) +
  geom_line(aes(y = predicted, color = "Forecast"), linewidth = 1) +
  geom_point(aes(y = predicted, color = "Forecast"), size = 2.5) +
  geom_line(aes(y = observed, color = "Observed"), linewidth = 1,
            linetype = "dashed") +
  geom_point(aes(y = observed, color = "Observed"), size = 2.5) +
  scale_x_continuous(breaks = seq_len(forecast_horizon),
                     labels = paste0("W", 201:207)) +
  scale_color_manual(values = c(Forecast = "blue", Observed = "red")) +
  labs(x = "Study week (Nov 1–Dec 19, 2025)",
       y = "Weekly dengue cases", color = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
print(final_plot)
ggsave(file.path(out_dir, "seven_week_forecast_with_PI.png"),
       final_plot, width = 10, height = 6, dpi = 300)
cat("\nAll outputs saved to: ", out_dir, "\n")
