# ============================================================
# Autocorrelation Diagnostics: (a) ACF, (b) PACF, (c) Ljung-Box
# ============================================================

library(dplyr)
library(ggplot2)
library(forecast)
library(patchwork)   # install.packages("patchwork") if not already installed

# -----------------------------
# Publication theme
# -----------------------------
theme_pub <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", linewidth = 0.5),
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# -----------------------------
# Load and order data
# -----------------------------
data1 <- read.csv("MyBDweeklyFinaldata.csv",stringsAsFactors = FALSE)
data1 <- data1[order(data1$Week), ]

max_lag <- 24

# -----------------------------
# (a) ACF panel
# -----------------------------
p1 <- ggAcf(data1$Dengue_cases, lag.max = max_lag) +
  labs(x = "Lag (weeks)", y = "Autocorrelation", title = "(a) ACF") +
  theme_pub()
print(p1)
# -----------------------------
# (b) PACF panel
# -----------------------------
p2 <- ggPacf(data1$Dengue_cases, lag.max = max_lag) +
  labs(x = "Lag (weeks)", y = "Partial autocorrelation", title = "(b) PACF") +
  theme_pub()
print(p2)
# -----------------------------
# (c) Ljung-Box statistic across lags 1-24
# -----------------------------
ljung_box_results <- lapply(1:max_lag, function(lag_value) {
  test <- Box.test(data1$Dengue_cases, lag = lag_value, type = "Ljung-Box")
  data.frame(
    Lag = lag_value,
    Statistic = as.numeric(test$statistic),
    p_value = test$p.value
  )
}) %>%
  bind_rows()

cat("\n--- Ljung-Box test results (lags 1-24) ---\n")
print(ljung_box_results)
write.csv(ljung_box_results, "ljung_box_results.csv", row.names = FALSE)

p3 <- ggplot(ljung_box_results, aes(x = Lag, y = Statistic)) +
  geom_line(colour = "#2C7FB8", linewidth = 1) +
  geom_point(colour = "#2C7FB8", size = 2.5) +
  scale_x_continuous(breaks = seq(2, max_lag, by = 2), limits = c(1, max_lag)) +
  labs(
    x = "Lag (weeks)",
    y = expression(paste("Ljung-Box ", chi^2, " statistic")),
    title = "(c) Ljung-Box Statistic"
  ) +
  theme_pub(base_size = 14)
print(p3)


ggsave("ACF.pdf",   plot = p1, width = 8, height = 5, units = "in")
ggsave("PACF.pdf",  plot = p2, width = 8, height = 5, units = "in")
ggsave("LjungBox.png", plot = p3, width = 7, height = 5, units = "in", dpi = 600, bg = "white")
