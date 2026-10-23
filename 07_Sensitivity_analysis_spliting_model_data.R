
# ============================================================
# CONSOLIDATED DENGUE SENSITIVITY ANALYSIS
#
# A. FIXED BASELINE LAGS, DIFFERENT TEST WINDOWS
# B. RESELECTED LAGS, DIFFERENT TEST WINDOWS

# ============================================================
rm(list = ls())
library(dplyr)
library(splines)
library(ggplot2)

# ============================================================
# 1. LOAD AND PREPARE DATA
# ============================================================

df <- read.csv("MyBDweeklyFinaldata.csv",  stringsAsFactors = FALSE)

df <- df %>%
  dplyr::arrange(Week) %>%
  dplyr::mutate(
    time_index = dplyr::row_number(),
    trend = time_index,
    sin1 = sin(2 * pi * time_index / 52),
    cos1 = cos(2 * pi * time_index / 52)
  )

required_cols <- c(
  "Week", "Dengue_cases",
  "Temperature", "Precipitation", "Humidity"
)

stopifnot(
  all(required_cols %in% names(df)),
  nrow(df) == 200,
  identical(as.numeric(df$Week), as.numeric(1:200)),
  !anyNA(df[, required_cols]),
  all(df$Dengue_cases >= 0)
)

# ============================================================
# 2. MODEL SPECIFICATIONS
# ============================================================

# A: Fixed baseline model for all partitions

fixed_specs <- data.frame(
  Split = c("84/16", "80/20", "70/30"),
  train_end = c(168, 160, 140),
  ar_order = c(8, 8, 8),
  temp_lag = c(9, 9, 9),
  precip_lag = c(7, 7, 7),
  humid_lag = c(8, 8, 8)
)

# B and C: Independently selected lag structures

selected_specs <- data.frame(
  Split = c("84/16", "80/20", "70/30"),
  train_end = c(168, 160, 140),
  ar_order = c(8, 8, 1),
  temp_lag = c(9, 15, 10),
  precip_lag = c(7, 7, 7),
  humid_lag = c(8, 5, 7)
)

# ============================================================
# 3. PERFORMANCE FUNCTIONS
# ============================================================

rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}

mae <- function(obs, pred) {
  mean(abs(obs - pred), na.rm = TRUE)
}

# ============================================================
# 4. CONSTRUCT LAGGED VARIABLES
# ============================================================

make_model_data <- function(
    data, ar_order, temp_lag,
    precip_lag, humid_lag
) {
  
  out <- data
  
  for (i in seq_len(ar_order)) {
    out[[paste0("lag_dengue_", i)]] <-
      log(dplyr::lag(out$Dengue_cases, i) + 1)
  }
  
  out$Temp_lag <- dplyr::lag(
    out$Temperature, temp_lag
  )
  
  out$Precip_lag <- dplyr::lag(
    out$Precipitation, precip_lag
  )
  
  out$Humd_lag <- dplyr::lag(
    out$Humidity, humid_lag
  )
  
  out
}

# ============================================================
# 5. NATURAL SPLINE FUNCTIONS
# ============================================================

build_ns_basis <- function(train_x) {
  
  train_x <- train_x[is.finite(train_x)]
  
  knots <- as.numeric(
    quantile(
      train_x,
      probs = c(0.25, 0.50, 0.75),
      na.rm = TRUE
    )
  )
  
  knots <- unique(knots)
  
  boundary_knots <- range(
    train_x,
    na.rm = TRUE
  )
  
  basis_train <- splines::ns(
    train_x,
    knots = knots,
    Boundary.knots = boundary_knots
  )
  
  list(
    train_basis = basis_train,
    knots = knots,
    boundary_knots = boundary_knots,
    predict_fn = function(new_x) {
      predict(basis_train, new_x)
    }
  )
}

add_spline_terms <- function(
    train_data, test_data = NULL
) {
  
  temp_basis <- build_ns_basis(train_data$Temp_lag)
  precip_basis <- build_ns_basis(train_data$Precip_lag)
  humid_basis <- build_ns_basis(train_data$Humd_lag)
  
  train_temp <- as.data.frame(
    temp_basis$train_basis
  )
  
  train_prec <- as.data.frame(
    precip_basis$train_basis
  )
  
  train_humd <- as.data.frame(
    humid_basis$train_basis
  )
  
  names(train_temp) <- paste0(
    "ns_temp_", seq_len(ncol(train_temp))
  )
  
  names(train_prec) <- paste0(
    "ns_precip_", seq_len(ncol(train_prec))
  )
  
  names(train_humd) <- paste0(
    "ns_humid_", seq_len(ncol(train_humd))
  )
  
  train_full <- dplyr::bind_cols(
    train_data,
    train_temp,
    train_prec,
    train_humd
  )
  
  spline_terms <- c(
    names(train_temp),
    names(train_prec),
    names(train_humd)
  )
  
  if (is.null(test_data)) {
    return(list(
      train_full = train_full,
      spline_terms = spline_terms,
      temp_basis = temp_basis,
      precip_basis = precip_basis,
      humid_basis = humid_basis
    ))
  }
  
  test_temp <- as.data.frame(
    temp_basis$predict_fn(test_data$Temp_lag)
  )
  
  test_prec <- as.data.frame(
    precip_basis$predict_fn(test_data$Precip_lag)
  )
  
  test_humd <- as.data.frame(
    humid_basis$predict_fn(test_data$Humd_lag)
  )
  
  names(test_temp) <- names(train_temp)
  names(test_prec) <- names(train_prec)
  names(test_humd) <- names(train_humd)
  
  test_full <- dplyr::bind_cols(
    test_data,
    test_temp,
    test_prec,
    test_humd
  )
  
  list(
    train_full = train_full,
    test_full = test_full,
    spline_terms = spline_terms,
    temp_basis = temp_basis,
    precip_basis = precip_basis,
    humid_basis = humid_basis
  )
}

# ============================================================
# 6. GENERAL MODEL EVALUATION FUNCTION
# ============================================================

evaluate_model <- function(
    spec,
    data,
    analysis_name,
    common_window = FALSE
) {
  
  split_name <- as.character(spec$Split)
  
  train_end <- as.integer(spec$train_end)
  ar_order <- as.integer(spec$ar_order)
  temp_lag <- as.integer(spec$temp_lag)
  precip_lag <- as.integer(spec$precip_lag)
  humid_lag <- as.integer(spec$humid_lag)
  
  max_lag <- max(
    ar_order, temp_lag,
    precip_lag, humid_lag
  )
  
  # ----------------------------------------------------------
  # Determine testing period
  # ----------------------------------------------------------
  
  if (common_window) {
    test_start <- 169
  } else {
    test_start <- train_end + 1
  }
  
  test_end <- 200
  
  # ----------------------------------------------------------
  # Construct training data
  # ----------------------------------------------------------
  
  train_raw <- data %>%
    dplyr::filter(time_index <= train_end)
  
  train_lagged <- make_model_data(
    train_raw,
    ar_order,
    temp_lag,
    precip_lag,
    humid_lag
  )
  
  # ----------------------------------------------------------
  # Construct test data with sufficient lag history
  #
  # Earlier observed counts are permitted as lag inputs,
  # but are not included in model fitting.
  # ----------------------------------------------------------
  
  test_raw <- data %>%
    dplyr::filter(
      time_index >= test_start - max_lag
    )
  
  test_lagged <- make_model_data(
    test_raw,
    ar_order,
    temp_lag,
    precip_lag,
    humid_lag
  )
  
  ar_terms <- paste0(
    "lag_dengue_", seq_len(ar_order)
  )
  
  needed_cols <- c(
    "Dengue_cases",
    ar_terms,
    "Temp_lag",
    "Precip_lag",
    "Humd_lag",
    "trend",
    "sin1",
    "cos1"
  )
  
  train_model <- train_lagged %>%
    dplyr::filter(
      dplyr::if_all(
        dplyr::all_of(needed_cols),
        ~ !is.na(.)
      )
    )
  
  test_model <- test_lagged %>%
    dplyr::filter(
      dplyr::if_all(
        dplyr::all_of(needed_cols),
        ~ !is.na(.)
      )
    ) %>%
    dplyr::filter(
      time_index >= test_start,
      time_index <= test_end
    )
  
  stopifnot(
    max(train_model$time_index) == train_end,
    min(test_model$time_index) == test_start,
    max(test_model$time_index) == test_end,
    nrow(test_model) == test_end - test_start + 1
  )
  
  # ----------------------------------------------------------
  # Training-only spline construction
  # ----------------------------------------------------------
  
  spline_obj <- add_spline_terms(
    train_model,
    test_model
  )
  
  train_full <- spline_obj$train_full
  test_full <- spline_obj$test_full
  
  spline_terms <- spline_obj$spline_terms
  
  # ----------------------------------------------------------
  # Fit original quasi-Poisson model
  # ----------------------------------------------------------
  
  model_formula <- stats::reformulate(
    c(
      ar_terms,
      spline_terms,
      "trend",
      "sin1",
      "cos1"
    ),
    response = "Dengue_cases"
  )
  
  fit <- stats::glm(
    model_formula,
    family = stats::quasipoisson(link = "log"),
    data = train_full
  )
  
  if (!isTRUE(fit$converged)) {
    stop(paste("Nonconvergence:", analysis_name, split_name))
  }
  
  # ----------------------------------------------------------
  # Original prediction procedure
  # ----------------------------------------------------------
  
  train_pred <- stats::predict(
    fit,
    newdata = train_full,
    type = "response"
  )
  
  test_pred <- stats::predict(
    fit,
    newdata = test_full,
    type = "response"
  )
  
  stopifnot(
    all(is.finite(train_pred)),
    all(is.finite(test_pred))
  )
  
  # ----------------------------------------------------------
  # Original error metrics
  # ----------------------------------------------------------
  
  train_sd <- stats::sd(
    train_full$Dengue_cases,
    na.rm = TRUE
  )
  
  stopifnot(is.finite(train_sd), train_sd > 0)
  
  train_rmse <- rmse(
    train_full$Dengue_cases,
    train_pred
  )
  
  train_mae <- mae(
    train_full$Dengue_cases,
    train_pred
  )
  
  test_rmse <- rmse(
    test_full$Dengue_cases,
    test_pred
  )
  
  test_mae <- mae(
    test_full$Dengue_cases,
    test_pred
  )
  
  # ----------------------------------------------------------
  # Model summary
  # ----------------------------------------------------------
  
  summary_row <- data.frame(
    Analysis = analysis_name,
    Split = split_name,
    Train_start = 1,
    Train_end = train_end,
    Test_start = test_start,
    Test_end = test_end,
    AR = ar_order,
    Temp_lag = temp_lag,
    Precip_lag = precip_lag,
    RH_lag = humid_lag,
    n_train = nrow(train_full),
    n_test = nrow(test_full),
    Train_RMSE = train_rmse,
    Train_MAE = train_mae,
    Train_SRMSE = train_rmse / train_sd,
    Train_SMAE = train_mae / train_sd,
    Test_RMSE = test_rmse,
    Test_MAE = test_mae,
    Test_SRMSE = test_rmse / train_sd,
    Test_SMAE = test_mae / train_sd
  )
  
  # ----------------------------------------------------------
  # Predictions for visualization
  # ----------------------------------------------------------
  
  prediction_row <- data.frame(
    Analysis = analysis_name,
    Split = split_name,
    Week = test_full$time_index,
    Observed = test_full$Dengue_cases,
    Predicted = as.numeric(test_pred)
  )
  
  list(
    summary = summary_row,
    predictions = prediction_row,
    fit = fit
  )
}

# ============================================================
# 7. RUN ANALYSIS A: FIXED BASELINE LAGS
# ============================================================

run_analysis <- function(
    specs,
    analysis_name,
    common_window = FALSE
) {
  
  outputs <- lapply(
    seq_len(nrow(specs)),
    function(i) {
      evaluate_model(
        spec = specs[i, ],
        data = df,
        analysis_name = analysis_name,
        common_window = common_window
      )
    }
  )
  
  list(
    summary = dplyr::bind_rows(
      lapply(outputs, function(x) x$summary)
    ),
    predictions = dplyr::bind_rows(
      lapply(outputs, function(x) x$predictions)
    ),
    fits = lapply(outputs, function(x) x$fit)
  )
}

analysis_A <- run_analysis(
  fixed_specs,
  analysis_name = "A: Fixed baseline",
  common_window = FALSE
)

# ============================================================
# 8. RUN ANALYSIS B: RESELECTED LAGS
# ============================================================

analysis_B <- run_analysis(
  selected_specs,
  analysis_name = "B: Reselected",
  common_window = FALSE
)


# ============================================================
# 10. COMBINE ALL RESULTS
# ============================================================

all_results <- dplyr::bind_rows(
  analysis_A$summary,
  analysis_B$summary,
)

# ============================================================
# 11. FORMAT TABLES
# ============================================================

format_table <- function(x) {
  
  x %>%
    dplyr::transmute(
      Split,
      AR,
      Temp_lag,
      Precip_lag,
      RH_lag,
      n_train,
      n_test,
      Train_RMSE = round(Train_RMSE, 2),
      Test_RMSE = round(Test_RMSE, 2),
      Test_MAE = round(Test_MAE, 2),
      Test_SRMSE = round(Test_SRMSE, 4),
      Test_SMAE = round(Test_SMAE, 4)
    )
}

table_A <- format_table(analysis_A$summary)
table_B <- format_table(analysis_B$summary)

cat("\n============================================\n")
cat("A: FIXED BASELINE STRUCTURE\n")
cat("============================================\n")
print(as.data.frame(table_A), row.names = FALSE)

cat("\n============================================\n")
cat("B: INDEPENDENTLY RESELECTED STRUCTURES\n")
cat("============================================\n")
print(as.data.frame(table_B), row.names = FALSE)


# ============================================================
# 12. COMBINED COMPARISON TABLE
# ============================================================

combined_table <- all_results %>%
  dplyr::transmute(
    Analysis,
    Split,
    Train_end,
    Test_start,
    Test_end,
    AR,
    Temp_lag,
    Precip_lag,
    RH_lag,
    n_train,
    n_test,
    Train_RMSE = round(Train_RMSE, 2),
    Test_RMSE = round(Test_RMSE, 2),
    Test_MAE = round(Test_MAE, 2),
    Test_SRMSE = round(Test_SRMSE, 4),
    Test_SMAE = round(Test_SMAE, 4)
  )

cat("\n============================================\n")
cat("COMBINED SENSITIVITY ANALYSIS TABLE\n")
cat("============================================\n")

print(
  as.data.frame(combined_table),
  row.names = FALSE
)

