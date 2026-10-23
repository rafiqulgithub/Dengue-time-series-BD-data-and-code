# ============================================================
# CROSS-LAGGED SCATTER PLOTS: TRAINING DATA ONLY
# Training dataset: Weeks 1-168
# Holdout dataset:  Weeks 169-200 
# ============================================================

library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load and order data
# ------------------------------------------------------------


data1 <- read.csv("MyBDweeklyFinaldata.csv",stringsAsFactors = FALSE) %>%
  arrange(Week)

# ------------------------------------------------------------
# 2. Restrict to training data
# ------------------------------------------------------------

train <- data1 %>%
  filter(Week >= 1 & Week <= 168)

stopifnot(
  nrow(data1) == 200,
  nrow(train) == 168,
  identical(as.numeric(train$Week), as.numeric(1:168)),
  !anyNA(train$Dengue_cases)
)

cat("Total observations:", nrow(data1), "\n")
cat("Training observations:", nrow(train), "\n")
cat("Training weeks:", min(train$Week), "-", max(train$Week), "\n")
cat("Holdout excluded: TRUE\n")

# ------------------------------------------------------------
# 3. Function for cross-lagged scatter plots
# ------------------------------------------------------------

make_lag_scatter_plot <- function(
    data,
    var_name,
    var_label,
    lag_range,
    nrow,
    ncol,
    x_breaks = NULL
) {
  
  # Construct lagged meteorological variables
  lagged_list <- lapply(lag_range, function(l) {
    
    data.frame(
      Dengue_cases = data$Dengue_cases,
      Value = if (l == 0) {
        data[[var_name]]
      } else {
        dplyr::lag(data[[var_name]], l)
      },
      Lag = l
    )
    
  })
  
  scatter_data <- bind_rows(lagged_list) %>%
    filter(
      !is.na(Value),
      !is.na(Dengue_cases)
    )
  
  # Panel labels
  lag_names <- ifelse(
    lag_range == 0,
    paste0(var_label, "(t)"),
    paste0(var_label, "(t-", lag_range, ")")
  )
  
  scatter_data$Lag_label <- factor(
    ifelse(
      scatter_data$Lag == 0,
      paste0(var_label, "(t)"),
      paste0(var_label, "(t-", scatter_data$Lag, ")")
    ),
    levels = lag_names
  )
  
  # Pearson correlations
  cor_labels <- scatter_data %>%
    group_by(Lag_label) %>%
    summarise(
      r = cor(
        Value,
        Dengue_cases,
        use = "complete.obs"
      ),
      .groups = "drop"
    ) %>%
    mutate(
      label = sprintf("%.2f", r)
    )
  
  # Construct plot
  p <- ggplot(
    scatter_data,
    aes(x = Value, y = Dengue_cases)
  ) +
    
    geom_point(
      shape = 1,
      size = 1.5,
      alpha = 1.5,
      colour = "black"
    ) +
    
    geom_smooth(
      method = "loess",
      se = FALSE,
      colour = "red",
      linewidth = 1.25
    ) +
    
    geom_label(
      data = cor_labels,
      aes(
        x = Inf,
        y = Inf,
        label = label
      ),
      hjust = 1.2,
      vjust = 1.5,
      inherit.aes = FALSE,
      size = 3.2,
      colour = "black",
      fill = "white",
      linewidth = 0.5,
      label.padding = grid::unit(0.15, "lines")
    ) +
    
    facet_wrap(
      ~ Lag_label,
      nrow = nrow,
      ncol = ncol
    ) +
    
    labs(
      x = paste0(var_label, "(t-lag)"),
      y = "Dengue cases(t)"
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      panel.grid.major = element_line(
        colour = "grey85",
        linewidth = 0.5
      ),
      
      panel.grid.minor = element_line(
        colour = "grey92",
        linewidth = 0.3
      ),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1.2
      ),
      
      strip.background = element_blank(),
      
      # Hide facet headers, as in your original figure.
      # Individual lag labels appear on the x-axis instead.
      strip.text = element_blank(),
      
      axis.title = element_text(
        face = "bold"
      ),
      
      axis.text = element_text(
        colour = "black"
      )
    )
  
  if (!is.null(x_breaks)) {
    p <- p +
      scale_x_continuous(
        breaks = x_breaks
      )
  }
  
  return(p)
}

# ============================================================
# 4. TEMPERATURE: Lags 0-15, 4 x 4
# ============================================================

p_temp <- make_lag_scatter_plot(
  data = train,
  var_name = "Temperature",
  var_label = "Temperature",
  lag_range = 0:15,
  nrow = 4,
  ncol = 4
)

print(p_temp)

ggsave(
  filename = "temperature_cross_lagged_training.png",
  plot = p_temp,
  width = 12,
  height = 10,
  dpi = 600,
  bg = "white"
)

# ============================================================
# 5. HUMIDITY: Lags 0-8, 3 x 3
# ============================================================

p_humid <- make_lag_scatter_plot(
  data = train,
  var_name = "Humidity",
  var_label = "Humidity",
  lag_range = 0:8,
  nrow = 3,
  ncol = 3
)

print(p_humid)

ggsave(
  filename = "humidity_cross_lagged_training.png",
  plot = p_humid,
  width = 12,
  height = 9,
  dpi = 600,
  bg = "white"
)

# ============================================================
# 6. PRECIPITATION: Lags 0-8, 3 x 3
# ============================================================

p_precip <- make_lag_scatter_plot(
  data = train,
  var_name = "Precipitation",
  var_label = "Precipitation",
  lag_range = 0:8,
  nrow = 3,
  ncol = 3
)

print(p_precip)

ggsave(
  filename = "precipitation_cross_lagged_training.png",
  plot = p_precip,
  width = 12,
  height = 9,
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 7. Print correlation tables
# ------------------------------------------------------------

get_lag_correlations <- function(data, var_name, lags) {
  
  bind_rows(lapply(lags, function(l) {
    
    x <- dplyr::lag(data[[var_name]], l)
    y <- data$Dengue_cases
    
    valid <- complete.cases(x, y)
    
    data.frame(
      Variable = var_name,
      Lag = l,
      Correlation = cor(x[valid], y[valid]),
      N = sum(valid)
    )
    
  }))
  
}

correlation_results <- bind_rows(
  get_lag_correlations(train, "Temperature", 0:15),
  get_lag_correlations(train, "Humidity", 0:8),
  get_lag_correlations(train, "Precipitation", 0:8)
)

print(correlation_results)

write.csv(
  correlation_results,
  "cross_lagged_correlations_training.csv",
  row.names = FALSE
)

cat("\nAll training-only cross-lagged scatter plots completed.\n")