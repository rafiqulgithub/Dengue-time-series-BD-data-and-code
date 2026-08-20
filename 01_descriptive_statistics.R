# Load data
df <-  read.csv("MyBDweeklyFinaldata.csv",stringsAsFactors = FALSE)

# ---- Basic structure ----
str(df)
dim(df)
summary(df)

# ---- Check for missing values ----
colSums(is.na(df))

# ---- Descriptive statistics function ----
library(dplyr)

descriptive_stats <- function(x) {
  c(
    N       = length(x),
    Mean    = mean(x, na.rm = TRUE),
    SD      = sd(x, na.rm = TRUE),
    Median  = median(x, na.rm = TRUE),
    Min     = min(x, na.rm = TRUE),
    Max     = max(x, na.rm = TRUE),
    Q1      = quantile(x, 0.25, na.rm = TRUE),
    Q3      = quantile(x, 0.75, na.rm = TRUE),
    IQR     = IQR(x, na.rm = TRUE),
    CV      = sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100,
    Skew    = e1071::skewness(x, na.rm = TRUE),
    Kurtosis = e1071::kurtosis(x, na.rm = TRUE)
  )
}

# Install e1071 if not already installed (for skewness/kurtosis)
if (!require(e1071)) install.packages("e1071")
library(e1071)

vars <- c("Dengue_cases", "Temperature", "Humidity", "Precipitation")
desc_table <- t(sapply(df[vars], descriptive_stats))
desc_table <- round(desc_table, 3)
print(desc_table)

# ---- Correlation matrix ----
cor_matrix <- cor(df[vars], use = "complete.obs", method = "pearson")
print(round(cor_matrix, 3))

# Optional: correlation with significance (p-values)
library(Hmisc)
cor_results <- rcorr(as.matrix(df[vars]))
print(cor_results$r)   # correlation coefficients
print(cor_results$P)   # p-values




# Time series plots (since data is weekly)
par(mfrow = c(2, 2))
for (v in vars) {
  plot(df$Week, df[[v]], type = "l", col = "darkred",
       main = paste(v, "over Weeks"), xlab = "Week", ylab = v)
}
par(mfrow = c(1, 1))

# ---- Correlation plot (visual) ----
if (!require(corrplot)) install.packages("corrplot")
library(corrplot)
corrplot(cor_matrix, method = "circle", type = "upper",
         addCoef.col = "black", tl.col = "black", tl.srt = 45)

# ---- Export summary table to CSV (optional) ----
write.csv(desc_table, "descriptive_stats_summary.csv")

