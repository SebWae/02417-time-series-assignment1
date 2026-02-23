# Libraries ----
library(tidyverse)

# Read data ----
#! Perhaps you need to set the working directory!?


setwd("/Users/krusand/Documents/GitHub/02417-time-series-assignment1/")
#setwd("C:/Users/sebas/OneDrive/Documents/MSc_HCAI/2_time_series/02417-time-series-assignment1")
D <- read.csv("DST_BIL54.csv")
str(D)

# See the help
?strftime
D$time <- as.POSIXct(paste0(D$time,"-01"), "%Y-%m-%d", tz="UTC")
D$time
class(D$time)

## Year to month for each of them
D$year <- 1900 + as.POSIXlt(D$time)$year + as.POSIXlt(D$time)$mon / 12

## Make the output variable a floating point (i.e.\ decimal number)
D$total <- as.numeric(D$total) / 1E6

## Divide intro train and test set
teststart <- as.POSIXct("2024-01-01", tz="UTC")
Dtrain <- D[D$time < teststart, ]
Dtest <- D[D$time >= teststart, ]


# 1. Plot data ----

ggplot(data=Dtrain, aes(x=time, y=total)) +
  geom_line() +
  xlab("Time [monthly]") +
  ylab("Total cars registered in Denmark at month [millions]") +
  ggtitle("The total number of cars registered at a given month in Denmark vs. time")

# 2. Linear Trend Model ----

