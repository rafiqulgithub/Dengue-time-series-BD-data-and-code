
library(dplyr)
library(ggplot2)

df <-   read.csv("MyBDweeklyFinaldata.csv",stringsAsFactors = FALSE)

traindf<-data.frame(df[1:168,])
dim(traindf)
# View structure
str(traindf)
str(testdf)
#View(testdf)

Dengue_show_yearly=traindf %>%dplyr::select(Dengue_cases)
Dengue_show_yearly=ts(Dengue_show_yearly,freq=52,start =c(2022,1))
#################################################################
plot(Dengue_show_yearly,
     xaxt = "n",  # suppress default x-axis
     xlab = "",
     ylab = "Dengue Cases",
     main = "Weekly Dengue Cases (2022–2025)")

# total weeks
weeks <- length(Dengue_show_yearly)

# week labels (2022w1, 2022w2, ...)
week_labels <- paste0(2022 + ((0:(weeks-1)) %/% 52),
                      "w", ((0:(weeks-1)) %% 52) + 1)

# get the actual time positions from the ts object
time_points <- time(Dengue_show_yearly)

# show every 4th week (adjust 'by' to control density)
axis(1,
     at = time_points[seq(1, weeks, by = 2)],
     labels = week_labels[seq(1, weeks, by = 2)],
     las = 2,
     cex.axis = 0.7)

#################################################################


Temp_show_yearly=traindf %>%dplyr::select(Temperature)
Temp_show_yearly=ts(Temp_show_yearly,freq=52,start =c(2022,1))

Humid_show_yearly=traindf %>%dplyr::select(Humidity)
Humid_show_yearly=ts(Humid_show_yearly,freq=52,start =c(2022,1))


Precip_show_yearly=traindf %>%dplyr::select(Precipitation)
Precip_show_yearly=ts(Precip_show_yearly,freq=52,start =c(2022,1))


Temperature <- ts(as.numeric(Temp_show_yearly[, 1]),
                  start = start(Temp_show_yearly),
                  frequency = frequency(Temp_show_yearly))

Precipitation <- ts(as.numeric(Precip_show_yearly[, 1]),
                    start = start(Precip_show_yearly),
                    frequency = frequency(Precip_show_yearly))

Humidity <- ts(as.numeric(Humid_show_yearly[, 1]),
               start = start(Humid_show_yearly),
               frequency = frequency(Humid_show_yearly))

Dengue_cases <- ts(as.numeric(Dengue_show_yearly[, 1]),
                   start = start(Dengue_show_yearly),
                   frequency = frequency(Dengue_show_yearly))

#install.packages('astsa')
library(astsa)

astsa::lag2.plot(Temperature,Dengue_cases,15)
astsa::lag2.plot(Precipitation,Dengue_cases,8)
astsa::lag2.plot(Humidity,Dengue_cases,8)