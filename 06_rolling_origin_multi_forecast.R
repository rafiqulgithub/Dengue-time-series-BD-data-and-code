# ================================================================
# DENGUE: 26-ORIGIN RECURSIVE FORECAST VALIDATION
# + NON-OVERLAPPING SENSITIVITY ANALYSIS
# + FIVE-PANEL FIGURE (3 columns x 2 rows; 8 historical weeks)
# Historical forecast-error intervals
# ================================================================
rm(list = ls())
library(dplyr)
library(splines)
library(ggplot2)
library(patchwork)

# 0. SETTINGS ------------------------------------------------------
out_dir <- "~/Downloads/dengue_26_origin_with_nonoverlap"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

calibration_origins <- 100:161
evaluation_origins  <- 168:193
nonoverlap_origins  <- c(168L, 175L, 182L, 189L)
figure_origins      <- c(168L, 172L, 179L, 186L, 193L)

H <- 7L
ar_order <- 8L
temp_lag <- 9L
precip_lag <- 7L
humid_lag <- 8L
hist_weeks <- 8L
ncol_fig <- 3L
common_y <- FALSE

# 1. READ AND VALIDATE DATA ----------------------------------------
df <-  read.csv("MyBDweeklyFinaldata.csv", stringsAsFactors = FALSE)
required <- c("Week", "Dengue_cases", "Temperature",
              "Precipitation", "Humidity")
stopifnot(all(required %in% names(df)))
df <- df %>%
  dplyr::arrange(Week) %>%
  dplyr::mutate(time_index = dplyr::row_number(),
                trend = time_index,
                sin1 = sin(2*pi*time_index/52),
                cos1 = cos(2*pi*time_index/52))

stopifnot(
  nrow(df) == 200L,
  all(stats::complete.cases(df[, required])),
  all(is.finite(df$Dengue_cases)),
  all(df$Dengue_cases >= 0),
  min(temp_lag, precip_lag, humid_lag) >= H,
  max(calibration_origins) + H <= 168L,
  max(evaluation_origins) + H <= nrow(df),
  all(nonoverlap_origins %in% evaluation_origins),
  all(figure_origins %in% evaluation_origins),
  all(diff(nonoverlap_origins) >= H)
)
cases_all <- df$Dengue_cases

# 2. FIT MODEL UP TO EACH ORIGIN -----------------------------------
make_model_data <- function(dat) {
  out <- dat
  for (j in seq_len(ar_order)) {
    out[[paste0("lag_dengue_", j)]] <-
      log1p(dplyr::lag(out$Dengue_cases, j))
  }
  out$Temp_lag   <- dplyr::lag(out$Temperature, temp_lag)
  out$Precip_lag <- dplyr::lag(out$Precipitation, precip_lag)
  out$Humd_lag   <- dplyr::lag(out$Humidity, humid_lag)
  out
}

make_spline <- function(x) {
  knots <- unique(as.numeric(stats::quantile(
    x, c(0.25, 0.50, 0.75), na.rm = TRUE
  )))
  basis <- splines::ns(
    x, knots = knots, Boundary.knots = range(x)
  )
  list(training = basis,
       predict = function(z) stats::predict(basis, z))
}

fit_upto <- function(origin, model = c("full", "ar_only")) {
  model <- match.arg(model)
  lagged <- make_model_data(df[seq_len(origin), , drop = FALSE])
  ar_terms <- paste0("lag_dengue_", seq_len(ar_order))
  needed <- c("Dengue_cases", ar_terms, "trend", "sin1", "cos1")
  if (model == "full") {
    needed <- c(needed, "Temp_lag", "Precip_lag", "Humd_lag")
  }
  train <- lagged[stats::complete.cases(
    lagged[, needed, drop = FALSE]), , drop = FALSE]
  if (nrow(train) < 50L) {
    stop("Too few complete training rows at origin ", origin)
  }

  bases <- NULL
  spline_terms <- character()
  if (model == "full") {
    bases <- list(temp   = make_spline(train$Temp_lag),
                  precip = make_spline(train$Precip_lag),
                  humid  = make_spline(train$Humd_lag))
    for (nm in names(bases)) {
      b <- as.data.frame(bases[[nm]]$training)
      names(b) <- paste0("ns_", nm, "_", seq_len(ncol(b)))
      train <- dplyr::bind_cols(train, b)
      spline_terms <- c(spline_terms, names(b))
    }
  }

  f <- stats::reformulate(
    c(ar_terms, spline_terms, "trend", "sin1", "cos1"),
    response = "Dengue_cases"
  )
  fit <- stats::glm(f, family = stats::quasipoisson(link = "log"),
                    data = train)
  if (anyNA(stats::coef(fit))) {
    stop("Rank-deficient fit at origin ", origin, " (", model, ")")
  }
  if (!isTRUE(fit$converged)) {
    stop("GLM did not converge at origin ", origin)
  }
  list(fit = fit, bases = bases)
}

# 3. SEVEN-WEEK RECURSIVE FORECAST -------------------------------
recursive_forecast <- function(origin, model = c("full", "ar_only")) {
  model <- match.arg(model)
  obj <- fit_upto(origin, model)
  history <- tail(df$Dengue_cases[seq_len(origin)], ar_order)
  predictions <- numeric(H)

  for (h in seq_len(H)) {
    target <- origin + h
    newrow <- data.frame(
      trend = target,
      sin1 = sin(2*pi*target/52),
      cos1 = cos(2*pi*target/52)
    )
    for (j in seq_len(ar_order)) {
      newrow[[paste0("lag_dengue_", j)]] <-
        log1p(history[ar_order + 1L - j])
    }
    if (model == "full") {
      weather_idx <- c(
        temp = target - temp_lag,
        precip = target - precip_lag,
        humid = target - humid_lag
      )
      if (any(weather_idx < 1L) || any(weather_idx > origin)) {
        stop("Future weather required at origin ", origin)
      }
      weather <- list(
        temp = df$Temperature[weather_idx["temp"]],
        precip = df$Precipitation[weather_idx["precip"]],
        humid = df$Humidity[weather_idx["humid"]]
      )
      for (nm in names(weather)) {
        b <- as.data.frame(obj$bases[[nm]]$predict(weather[[nm]]))
        names(b) <- paste0("ns_", nm, "_", seq_len(ncol(b)))
        newrow <- dplyr::bind_cols(newrow, b)
      }
    }
    mu <- as.numeric(stats::predict(
      obj$fit, newdata = newrow, type = "response"
    ))
    if (length(mu) != 1L || !is.finite(mu) || mu < 0) {
      stop("Invalid forecast at origin ", origin, ", horizon ", h)
    }
    predictions[h] <- mu
    # For h > 1, use previously predicted means, not future observations.
    history <- c(history[-1], mu)
  }

  data.frame(
    model = model,
    origin = origin,
    horizon = seq_len(H),
    target_week = origin + seq_len(H),
    observed = cases_all[origin + seq_len(H)],
    predicted = predictions
  )
}

run_origins <- function(origins, phase,
                        models = c("full", "ar_only")) {
  dplyr::bind_rows(lapply(models, function(m) {
    dplyr::bind_rows(lapply(origins, function(o) {
      cat(phase, "|", m, "| origin", o, "\n")
      recursive_forecast(o, m)
    }))
  }))
}

# 4. CALIBRATION: DEVELOPMENT PERIOD ONLY -------------------------
calibration <- run_origins(calibration_origins, "Calibration") %>%
  dplyr::mutate(error = observed - predicted)

# Fixed historical signed-error quantiles for each model/horizon.
# These intervals are approximate; they have no formal 95% guarantee.
error_quantiles <- calibration %>%
  dplyr::group_by(model, horizon) %>%
  dplyr::summarise(
    n_calibration = dplyr::n(),
    error_q025 = unname(stats::quantile(error, 0.025, type = 1)),
    error_q975 = unname(stats::quantile(error, 0.975, type = 1)),
    calibration_MAE = mean(abs(error)),
    .groups = "drop"
  )

# 5. PRIMARY: 26 ORIGINS x 7 HORIZONS ------------------------------
evaluation <- run_origins(evaluation_origins, "Evaluation") %>%
  dplyr::left_join(error_quantiles, by = c("model", "horizon")) %>%
  dplyr::mutate(
    lower95 = pmax(0, predicted + error_q025),
    upper95 = pmax(0, predicted + error_q975),
    error = observed - predicted,
    absolute_error = abs(error),
    squared_error = error^2,
    ape = ifelse(observed > 0, 100 * absolute_error / observed, NA_real_),
    covered95 = observed >= lower95 & observed <= upper95,
    interval_width = upper95 - lower95
  )

stopifnot(
  nrow(evaluation) == 2L * length(evaluation_origins) * H,
  all(evaluation$lower95 <= evaluation$upper95),
  all(!is.na(evaluation$observed))
)

summarise_performance <- function(dat) {
  dat %>% dplyr::summarise(
    n = dplyr::n(),
    RMSE = sqrt(mean(squared_error)),
    MAE = mean(absolute_error),
    MAPE = mean(ape, na.rm = TRUE),
    Mean_Error = mean(error),
    Covered = sum(covered95),
    Coverage_percent = 100 * mean(covered95),
    Mean_PI_width = mean(interval_width),
    .groups = "drop"
  )
}

by_origin <- evaluation %>%
  dplyr::group_by(model, origin) %>% summarise_performance()
by_horizon <- evaluation %>%
  dplyr::group_by(model, horizon) %>% summarise_performance()
overall <- evaluation %>%
  dplyr::group_by(model) %>% summarise_performance()

# 6. SENSITIVITY: FOUR NON-OVERLAPPING SEVEN-WEEK WINDOWS ----------
# Origin 168 -> weeks 169:175; 175 -> 176:182;
# 182 -> 183:189; 189 -> 190:196.
# Reuse identical forecasts and intervals from the primary analysis.
sensitivity <- evaluation %>%
  dplyr::filter(origin %in% nonoverlap_origins)
sensitivity_by_origin <- sensitivity %>%
  dplyr::group_by(model, origin) %>% summarise_performance()
sensitivity_by_horizon <- sensitivity %>%
  dplyr::group_by(model, horizon) %>% summarise_performance()
sensitivity_overall <- sensitivity %>%
  dplyr::group_by(model) %>% summarise_performance()

stopifnot(
  nrow(sensitivity) == 2L * length(nonoverlap_origins) * H,
  dplyr::n_distinct(sensitivity$target_week) == 4L * H
)

# Explicit dplyr namespace avoids select() masking errors.
comparison <- dplyr::bind_rows(
  overall %>% dplyr::mutate(analysis = "Primary: 26 rolling origins"),
  sensitivity_overall %>%
    dplyr::mutate(analysis = "Sensitivity: 4 disjoint origins")
) %>%
  dplyr::select(analysis, dplyr::everything())

cat("\n--- CALIBRATION ERROR QUANTILES ---\n")
print(error_quantiles, n = Inf)
cat("\n--- PRIMARY: OVERALL ---\n")
print(overall, n = Inf)
cat("\n--- PRIMARY: BY HORIZON ---\n")
print(by_horizon, n = Inf)
cat("\n--- SENSITIVITY: OVERALL ---\n")
print(sensitivity_overall, n = Inf)
cat("\n--- SENSITIVITY: BY HORIZON ---\n")
print(sensitivity_by_horizon, n = Inf)
cat("\n--- COMPARISON ---\n")
print(comparison, n = Inf)

# 7. SAVE CSV TABLES ----------------------------------------------
write_table <- function(x, name) {
  utils::write.csv(x, file.path(out_dir, name), row.names = FALSE)
}
write_table(calibration, "calibration_recursive_errors.csv")
write_table(error_quantiles, "calibration_error_quantiles.csv")
write_table(evaluation, "primary_26_origin_forecasts.csv")
write_table(by_origin, "primary_performance_by_origin.csv")
write_table(by_horizon, "primary_performance_by_horizon.csv")
write_table(overall, "primary_performance_overall.csv")
write_table(sensitivity, "sensitivity_4_origin_forecasts.csv")
write_table(sensitivity_by_origin, "sensitivity_performance_by_origin.csv")
write_table(sensitivity_by_horizon, "sensitivity_performance_by_horizon.csv")
write_table(sensitivity_overall, "sensitivity_performance_overall.csv")
write_table(comparison, "primary_vs_sensitivity_comparison.csv")

# 8. FIGURE SETTINGS ----------------------------------------------
fig_font <- "sans"
fig_width_mm <- 220
col_obs <- "black"
col_fc <- "#0072B2"
palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")

fig_theme <- ggplot2::theme_classic(
  base_size = 10, base_family = fig_font
) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 9, hjust = 0),
    plot.subtitle = ggplot2::element_text(
      size = 7.5, colour = "grey30",
      margin = ggplot2::margin(b = 3)
    ),
    axis.title = ggplot2::element_text(size = 8.5),
    axis.text = ggplot2::element_text(size = 7.5, colour = "black"),
    axis.line = ggplot2::element_line(linewidth = 0.4),
    axis.ticks = ggplot2::element_line(linewidth = 0.4),
    panel.grid.major.y = ggplot2::element_line(
      colour = "grey90", linewidth = 0.3
    ),
    legend.text = ggplot2::element_text(size = 8),
    legend.title = ggplot2::element_text(size = 8, face = "bold"),
    legend.key.width = grid::unit(1, "cm"),
    legend.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(4, 8, 4, 4)
  )

# Save PDF/PNG; TIFF is optional because Quartz on macOS may not
# support LZW compression. Each format is saved independently.
save_fig <- function(p, stem, w_mm, h_mm) {
  files <- c(pdf = ".pdf", png = ".png", tiff = ".tiff")
  for (fmt in names(files)) {
    path <- file.path(out_dir, paste0(stem, files[[fmt]]))
    tryCatch({
      if (fmt == "pdf") {
        ggplot2::ggsave(path, p, width = w_mm, height = h_mm,
                        units = "mm", device = grDevices::pdf)
      } else if (fmt == "png") {
        ggplot2::ggsave(path, p, width = w_mm, height = h_mm,
                        units = "mm", dpi = 600, bg = "white")
      } else {
        ggplot2::ggsave(path, p, width = w_mm, height = h_mm,
                        units = "mm", dpi = 600, bg = "white",
                        compression = "none")
      }
      message("Saved: ", path)
    }, error = function(e) {
      message("Could not save ", fmt, " figure: ", conditionMessage(e))
    })
  }
}

# 9. FIGURE 1: FIVE PANELS, 3 IN FIRST ROW, 2 IN SECOND -----------
full_eval <- evaluation %>% dplyr::filter(model == "full")
origins <- figure_origins
n_panels <- length(origins)
common_ymax <- max(
  c(full_eval$upper95, full_eval$observed, full_eval$predicted),
  na.rm = TRUE
) * 1.05

make_panel <- function(i) {
  o <- origins[i]
  res <- full_eval %>% dplyr::filter(origin == o)
  wk <- seq.int(o - hist_weeks + 1L, o + H)
  obs <- data.frame(time_index = wk, observed = cases_all[wk])
  fc_line <- dplyr::bind_rows(
    data.frame(time_index = o, predicted = cases_all[o]),
    res %>% dplyr::transmute(time_index = target_week, predicted)
  )
  subtitle <- sprintf(
    "%d/%d observed weeks within approximate 95%% PI",
    sum(res$covered95), nrow(res)
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(
      xintercept = o, linetype = "dotted",
      colour = "grey40", linewidth = 0.4
    ) +
    ggplot2::geom_ribbon(
      data = res,
      ggplot2::aes(
        x = target_week, ymin = lower95, ymax = upper95,
        fill = "95% prediction interval"
      ),
      alpha = 0.18
    ) +
    ggplot2::geom_line(
      data = obs,
      ggplot2::aes(x = time_index, y = observed,
                   colour = "Observed", linetype = "Observed"),
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data = obs,
      ggplot2::aes(x = time_index, y = observed,
                   colour = "Observed", shape = "Observed"),
      size = 1.5
    ) +
    ggplot2::geom_line(
      data = fc_line,
      ggplot2::aes(x = time_index, y = predicted,
                   colour = "Forecast", linetype = "Forecast"),
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      data = res,
      ggplot2::aes(x = target_week, y = predicted,
                   colour = "Forecast", shape = "Forecast"),
      size = 1.7
    ) +
    ggplot2::annotate(
      "text", x = o, y = Inf, label = "Origin",
      hjust = 1.15, vjust = 1.6, size = 2.4, colour = "grey35"
    ) +
    ggplot2::scale_colour_manual(
      NULL, values = c(Observed = col_obs, Forecast = col_fc),
      breaks = c("Observed", "Forecast")
    ) +
    ggplot2::scale_linetype_manual(
      NULL, values = c(Observed = "solid", Forecast = "longdash"),
      breaks = c("Observed", "Forecast")
    ) +
    ggplot2::scale_shape_manual(
      NULL, values = c(Observed = 16, Forecast = 17),
      breaks = c("Observed", "Forecast")
    ) +
    ggplot2::scale_fill_manual(
      NULL, values = c("95% prediction interval" = col_fc)
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_pretty(n = 5),
      expand = ggplot2::expansion(mult = c(0.02, 0.03))
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0, 0.06))
    ) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(
      title = sprintf("(%s) Origin: week %d (forecast %d-%d)",
                      letters[i], o, o+1L, o+H),
      subtitle = subtitle,
      x = "Week index",
      y = "Weekly dengue cases"
    ) +
    fig_theme

  if (common_y) {
    p <- p + ggplot2::coord_cartesian(ylim = c(0, common_ymax))
  }
  p
}

panels <- lapply(seq_len(n_panels), make_panel)
# The sixth cell is reserved for the collected legend.
panels <- c(panels, list(patchwork::guide_area()))
fig1 <- patchwork::wrap_plots(panels, ncol = 3, nrow = 2) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
print(fig1)
save_fig(fig1, "Figure1_five_origins_3x2", fig_width_mm, 145)

