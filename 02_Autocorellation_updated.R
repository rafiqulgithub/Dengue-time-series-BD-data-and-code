# ============================================================
# Autocorrelation Diagnostics Using TRAINING DATA ONLY
# (a) ACF, (b) PACF, (c) Ljung-Box
# Training: Weeks 1-168
# Holdout:  Weeks 169-200 (NOT USED)
# ============================================================

library(dplyr)
library(ggplot2)
library(forecast)
library(patchwork)

# ------------------------------------------------------------
# 1. Publication theme
# ------------------------------------------------------------

theme_pub <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(
        colour = "black", linewidth = 1.5
      ),
      strip.background = element_rect(
        fill = "grey90", colour = "black"
      ),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# ------------------------------------------------------------
# 2. Load and order data
# ------------------------------------------------------------

data1 <- read.csv("MyBDweeklyFinaldata.csv",  stringsAsFactors = FALSE)

data1 <- data1 %>%
  arrange(Week)

# ------------------------------------------------------------
# 3. Split training and holdout data
# ------------------------------------------------------------

train <- data1 %>%
  filter(Week >= 1 & Week <= 168)

holdout <- data1 %>%
  filter(Week >= 169 & Week <= 200)

# Verify split
stopifnot(
  nrow(data1) == 200,
  nrow(train) == 168,
  nrow(holdout) == 32,
  identical(as.numeric(train$Week), as.numeric(1:168)),
  identical(as.numeric(holdout$Week), as.numeric(169:200)),
  !anyNA(train$Dengue_cases)
)

cat("\n========================================\n")
cat("DATA SPLIT VERIFICATION\n")
cat("========================================\n")

cat("Total observations:", nrow(data1), "\n")
cat("Training observations:", nrow(train), "\n")
cat("Training weeks:",
    min(train$Week), "-", max(train$Week), "\n")

cat("Holdout observations:", nrow(holdout), "\n")
cat("Holdout weeks:",
    min(holdout$Week), "-", max(holdout$Week), "\n")

cat("Holdout used in diagnostics: FALSE\n")

# ------------------------------------------------------------
# 4. Training dengue series
# ------------------------------------------------------------

y_train <- train$Dengue_cases

max_lag <- 24

# ------------------------------------------------------------
# 5. (a) ACF
# ------------------------------------------------------------

p1 <- ggAcf(
  y_train,
  lag.max = max_lag
) +
  labs(
    x = "Lag (weeks)",
    y = "Autocorrelation",
    title = "Autocorrelation Function(ACF)"
  ) +
  theme_pub()

print(p1)

# ------------------------------------------------------------
# 6. (b) PACF
# ------------------------------------------------------------

p2 <- ggPacf(
  y_train,
  lag.max = max_lag
) +
  labs(
    x = "Lag (weeks)",
    y = "Partial autocorrelation",
    title = "Partial autocorrelation Function(PACF)"
  ) +
  theme_pub()

print(p2)

# ------------------------------------------------------------
# 7. (c) Ljung-Box test across lags 1-24
# ------------------------------------------------------------

ljung_box_results <- lapply(
  1:max_lag,
  function(lag_value) {
    
    test <- Box.test(
      y_train,
      lag = lag_value,
      type = "Ljung-Box"
    )
    
    data.frame(
      Lag = lag_value,
      Statistic = as.numeric(test$statistic),
      p_value = test$p.value
    )
  }
) %>%
  bind_rows()

cat("\n========================================\n")
cat("LJUNG-BOX RESULTS: TRAINING DATA ONLY\n")
cat("========================================\n")

print(ljung_box_results)

write.csv(
  ljung_box_results,
  "ljung_box_results_training.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Ljung-Box statistic plot
# ------------------------------------------------------------

p3 <- ggplot(
  ljung_box_results,
  aes(x = Lag, y = Statistic)
) +
  geom_line(
    colour = "#2C7fb8",
    linewidth = 2.
  ) +
  geom_point(
    colour = "#2C7FB8",
    size = 2.6
  ) +
  scale_x_continuous(
    breaks = seq(2, max_lag, by = 2),
    limits = c(1, max_lag)
  ) +
  labs(
    x = "Lag (weeks)",
    y = expression(
      paste("Ljung-Box ", chi^2, " statistic")
    ),
    title = "Ljung-Box Statistic"
  ) +
  theme_pub(base_size = 14)

print(p3)

# ------------------------------------------------------------
# 9. Save publication figures
# ------------------------------------------------------------

ggsave(
  "ACF_training.png",
  plot = p1,
  width = 8,
  height = 5,
  units = "in"
)

ggsave(
  "PACF_training.png",
  plot = p2,
  width = 8,
  height = 5,
  units = "in"
)

ggsave(
  "LjungBox_training.png",
  plot = p3,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

cat("\nAll training-only diagnostics completed.\n")
