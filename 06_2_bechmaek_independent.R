
# ==============================================================================

library(dplyr)
library(splines)
library(MASS)
library(forecast)
library(knitr)

# ------------------------------------------------------------------------------
# 1. Load Data & Create Predictors
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

# Add AR(8) log-transformed dengue lags
for (i in 1:8) {
  df[[paste0("lag_dengue_", i)]] <- log(dplyr::lag(df$Dengue_cases, i) + 1)
}

# Split into Train / Test Sets
train_end_week <- 168

# Regression Training Set: Explicit complete-case filter for lagged features
train_data <- df %>% 
  filter(time_index <= train_end_week) %>% 
  filter(
    !is.na(lag_dengue_8),
    !is.na(Temp_lag),
    !is.na(Precip_lag),
    !is.na(Humd_lag)
  )

# Regression Holdout Test Set: Explicit complete-case filter for lagged features
test_data <- df %>% 
  filter(time_index > train_end_week) %>% 
  filter(
    !is.na(lag_dengue_8),
    !is.na(Temp_lag),
    !is.na(Precip_lag),
    !is.na(Humd_lag)
  )

# Full Development Set for Univariate SARIMA (Weeks 1-168)
sarima_train_data <- df %>%
  filter(time_index <= train_end_week)

# ------------------------------------------------------------------------------
# 2. Build Natural Spline Bases (Training Knots at 25th, 50th, 75th Percentiles)
# ------------------------------------------------------------------------------
build_ns_basis <- function(train_x) {
  x_clean <- train_x[is.finite(train_x)]
  knots <- as.numeric(quantile(x_clean, probs = c(0.25, 0.50, 0.75), na.rm = TRUE))
  knots <- unique(knots)
  boundary_knots <- range(x_clean, na.rm = TRUE)
  
  basis_train <- ns(train_x, knots = knots, Boundary.knots = boundary_knots)
  
  list(
    train_basis = basis_train,
    predict_fn  = function(new_x) predict(basis_train, new_x)
  )
}

temp_basis   <- build_ns_basis(train_data$Temp_lag)
precip_basis <- build_ns_basis(train_data$Precip_lag)
humid_basis  <- build_ns_basis(train_data$Humd_lag)

# Extract train and test spline matrices
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

colnames(train_splines) <- paste0("ns_col_", 1:ncol(train_splines))
colnames(test_splines)  <- colnames(train_splines)

train_full <- bind_cols(train_data, train_splines)
test_full  <- bind_cols(test_data, test_splines)

# Model Formulas
ar_terms <- paste0("lag_dengue_", 1:8, collapse = " + ")
spline_terms <- paste(colnames(train_splines), collapse = " + ")

formula_proposed <- as.formula(paste("Dengue_cases ~", ar_terms, "+", spline_terms, "+ trend + sin1 + cos1"))
formula_ar_only  <- as.formula(paste("Dengue_cases ~", ar_terms, "+ trend + sin1 + cos1"))

# ------------------------------------------------------------------------------
# 3. Model Fitting & 1-Step-Ahead Out-of-Sample Predictions
# ------------------------------------------------------------------------------

# Model 1: Proposed quasi-Poisson
fit_proposed <- glm(formula_proposed, family = quasipoisson(link = "log"), data = train_full)
pred_proposed <- predict(fit_proposed, newdata = test_full, type = "response")

# Model 2: Negative Binomial Regression
fit_nb <- tryCatch(
  MASS::glm.nb(formula_proposed, data = train_full),
  error = function(e) {
    warning("Negative Binomial model failed to converge: ", conditionMessage(e))
    NULL
  }
)
pred_nb <- if (!is.null(fit_nb)) predict(fit_nb, newdata = test_full, type = "response") else rep(NA_real_, nrow(test_full))

# Model 3: AR(8) quasi-Poisson (No Meteorology)
fit_ar_only <- glm(formula_ar_only, family = quasipoisson(link = "log"), data = train_full)
pred_ar_only <- predict(fit_ar_only, newdata = test_full, type = "response")

# Model 4: Sequential 1-Step-Ahead SARIMA
# Step A: Fit optimal SARIMA specification on the complete 168-week development set
ts_train <- ts(sarima_train_data$Dengue_cases, frequency = 52)
fit_sarima_train <- forecast::auto.arima(
  ts_train,
  seasonal = TRUE,
  stepwise = FALSE,
  approximation = FALSE
)


history_vec <- sarima_train_data$Dengue_cases
pred_sarima <- numeric(nrow(test_data))

for (i in seq_len(nrow(test_data))) {
  history_ts <- ts(history_vec, frequency = 52)
  
  fit_iter <- forecast::Arima(
    history_ts, 
    model = fit_sarima_train
  )
  
  pred_sarima[i] <- as.numeric(
    forecast::forecast(fit_iter, h = 1)$mean[1]
  )
  
  history_vec <- c(history_vec, test_data$Dengue_cases[i])
}

# Model 5: Seasonal Naïve Baseline (Y_t = Y_{t-52})
pred_snaive <- sapply(test_full$time_index, function(ti) {
  val <- df$Dengue_cases[df$time_index == (ti - 52)]
  if (length(val) == 1) val else NA_real_
})

# ------------------------------------------------------------------------------
# 4. Compute Performance Metrics
# ------------------------------------------------------------------------------
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae  <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)

train_sd <- sd(train_data$Dengue_cases, na.rm = TRUE)
obs_test <- test_full$Dengue_cases

calc_metrics <- function(obs, pred) {
  data.frame(
    RMSE  = round(rmse(obs, pred), 2),
    MAE   = round(mae(obs, pred), 2),
    SRMSE = round(rmse(obs, pred) / train_sd, 3),
    SMAE  = round(mae(obs, pred) /train_sd, 3)
  )
}

results <- bind_rows(
  cbind(Model = "Proposed quasi-Poisson (AR + Met Splines)", calc_metrics(obs_test, pred_proposed)),
  cbind(Model = "Negative Binomial (AR + Met Splines)", calc_metrics(obs_test, pred_nb)),
  cbind(Model = "AR(8) quasi-Poisson (No Meteorology)", calc_metrics(obs_test, pred_ar_only)),
  cbind(Model = "Sequential SARIMA", calc_metrics(obs_test, pred_sarima)),
  cbind(Model = "Seasonal Naïve Baseline", calc_metrics(obs_test, pred_snaive))
) %>% arrange(RMSE)

# Output Summary Table
knitr::kable(results, caption = "Out-of-Sample Performance Comparison on Holdout Test Set (One-Step-Ahead Evaluation)")
arimaorder(fit_sarima_train)