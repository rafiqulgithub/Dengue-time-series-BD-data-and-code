

required_packages <- c("dplyr", "ggplot2", "patchwork")
new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]
if (length(new_packages) > 0) install.packages(new_packages)

library(dplyr)
library(ggplot2)
library(patchwork)

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
# Output directory -- change to your local folder
# -----------------------------
out_dir <- "~/Downloads/ccf_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load data
# -----------------------------

data1 <-  read.csv("MyBDweeklyFinaldata.csv",stringsAsFactors = FALSE)
  arrange(Week) %>%
  mutate(time_index = row_number())

max_lag <- 24
n <- nrow(data1)

# Approximate 95% CCF significance bound
sig_bound <- 1.96 / sqrt(n)

# -----------------------------
# CCF helper function
# Positive lag = exposure leads dengue cases by k weeks
# (sign-corrected relative to R's native ccf() output; see note above)
# -----------------------------
run_ccf_lagwindow <- function(exposure, exposure_label, outcome, max_lag = 24) {
  
  cc <- ccf(
    x = exposure,
    y = outcome,
    lag.max = max_lag,
    plot = FALSE,
    na.action = na.pass
  )
  
  cc_df <- data.frame(
    lag = -cc$lag[, 1, 1],      # SIGN FLIP -- see note at top of script
    correlation = cc$acf[, 1, 1]
  ) %>%
    filter(lag >= 0, lag <= max_lag) %>%
    arrange(lag) %>%
    mutate(
      Variable = exposure_label,
      sig_upper = sig_bound,
      sig_lower = -sig_bound,
      significant = abs(correlation) > sig_bound
    )
  
  return(cc_df)
}

# -----------------------------
# Run CCF for each meteorological variable
# -----------------------------
ccf_temp <- run_ccf_lagwindow(
  exposure = data1$Temperature, exposure_label = "Temperature",
  outcome = data1$Dengue_cases, max_lag = max_lag
)
ccf_precip <- run_ccf_lagwindow(
  exposure = data1$Precipitation, exposure_label = "Precipitation",
  outcome = data1$Dengue_cases, max_lag = max_lag
)
ccf_humid <- run_ccf_lagwindow(
  exposure = data1$Humidity, exposure_label = "Humidity",
  outcome = data1$Dengue_cases, max_lag = max_lag
)

ccf_all <- bind_rows(ccf_temp, ccf_precip, ccf_humid)

print(ccf_all)
write.csv(ccf_all, file.path(out_dir, "ccf_lag_correlations.csv"), row.names = FALSE)

# -----------------------------
# Candidate window extraction: contiguous run of lags within
# 90% of the peak correlation (magnitude-based, not a
# significance threshold -- see earlier discussion on why the
# significance bound is uninformative at this sample size)
# -----------------------------
get_contiguous_window <- function(sig_lags, peak_lag) {
  if (length(sig_lags) == 0) return(c(NA, NA))
  sig_lags <- sort(sig_lags)
  runs <- split(sig_lags, cumsum(c(1, diff(sig_lags) != 1)))
  run_with_peak <- Filter(function(r) peak_lag %in% r, runs)
  if (length(run_with_peak) == 0) {
    dists <- sapply(runs, function(r) min(abs(r - peak_lag)))
    run_with_peak <- list(runs[[which.min(dists)]])
  }
  r <- run_with_peak[[1]]
  c(min(r), max(r))
}

ccf_summary <- ccf_all %>%
  filter(lag > 0) %>%   # exclude lag 0 -- no plausible same-week delay
  group_by(Variable) %>%
  summarise(
    strongest_lag = lag[which.max(correlation)],
    max_correlation = round(max(correlation), 3),
    candidate_window = {
      peak_val <- max(correlation)
      peak_lag_here <- lag[which.max(correlation)]
      prop_of_peak <- 0.90
      near_peak_lags <- lag[correlation >= prop_of_peak * peak_val]
      win <- get_contiguous_window(near_peak_lags, peak_lag_here)
      if (all(is.na(win))) "None" else paste0(win[1], "-", win[2], " weeks")
    },
    .groups = "drop"
  )

print(ccf_summary)
write.csv(ccf_summary, file.path(out_dir, "ccf_candidate_windows_summary.csv"), row.names = FALSE)

# -----------------------------
# Plot CCF curves
# -----------------------------
ccf_plot <- ggplot(ccf_all, aes(x = lag, y = correlation)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  geom_hline(yintercept = sig_bound, linetype = "dashed", colour = "#2C7FB8") +
  geom_hline(yintercept = -sig_bound, linetype = "dashed", colour = "#2C7FB8") +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(~ Variable, scales = "free_y", nrow = 1) +
  scale_x_continuous(breaks = seq(0, max_lag, by = 2), limits = c(0, max_lag)) +
  labs(x = "Lag (weeks)", y = "Cross-correlation coefficient") +
  theme_pub(base_size = 14)

print(ccf_plot)

ggsave(
  filename = file.path(out_dir, "CCF_lag_correlation_curves.png"),
  plot = ccf_plot, width = 10, height = 4, units = "in", dpi = 600, bg = "white"
)

cat("\nCCF exploratory analysis completed (sign-corrected).\n")
print(ccf_summary)

