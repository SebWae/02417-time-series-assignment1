# Libraries ----
library(tidyverse)

# Read data ----
#! Perhaps you need to set the working directory!?

# setwd("/Users/krusand/OneDrive/DTU/2. Semester/02417 - Time Series Analysis/Assignments/Assignment1/")
setwd("C:/Users/sebas/OneDrive/Documents/MSc_HCAI/2_time_series/02417-time-series-assignment1")
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


# 3. OLS - global linear trend model ----
## 3.1 Parameter estimation ----
# constructing the design matrix
X <- cbind(1, Dtrain$year)
print(X)

# constructing the y vector
y <- cbind(Dtrain$total)
print(y)

# estimating the parameters by solving the "normal equations":
theta_OLS <- solve(t(X)%*%X)%*%t(X)%*%y
print(theta_OLS)

# saving the parameter estimates as variables
theta_OLS_0 <- theta_OLS[1]
theta_OLS_1 <- theta_OLS[2]


## 3.2 Estimated standard errors and mean ----
# computing y_hat values (prediction)
yhat_OLS <- X%*%theta_OLS
Dtrain$yhat <- yhat_OLS

# computing residuals
e_OLS <- y - yhat_OLS

# calculating the sum of squared residuals
RSS_OLS <- t(e_OLS)%*%e_OLS

# we have 2 parameters theta_OLS_0 and theta_OLS_1
p <- 2

# n is the number of observations
n <- length(Dtrain$time)

# calculate sigma^2:
sigma2_OLS <- as.numeric(RSS_OLS/(n - p))

# calculate variance-covariance matrix of parameters:
V_ols <- sigma2_OLS * solve(t(X) %*% X)
print(V_ols)

# the variances of the parameters are the values in the diagonal:
diag(V_ols)
# and the standard errors are given by:
sqrt(diag(V_ols))


# plot estimated mean
ggplot(data=Dtrain, aes(x=time, y=yhat)) +
  geom_line() +
  xlab("Time [monthly]") +
  ylab("Number of cars registered [millions]") +
  ggtitle("Estimated number of cars registered in Denmark") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_point(aes(x=time, y=total))


## 3.3 Forecast with prediction interval ----
# design matrix for test set
X_test <- cbind(1, Dtest$year)

# predictions for the test set
yhat_test <- X_test%*%theta_OLS
Dtest$yhat <- yhat_test

# compute prediction variance-covariance matrix:
Vmatrix_pred <- sigma2_OLS*(1+(X_test%*%solve(t(X)%*%X))%*%t(X_test))
# the variances of individual predictions are in the diagonal of the matrix above
print(diag(Vmatrix_pred))

# compute "prediction intervals" (95% interval)
diff <- qt(0.975, df=n-p)*sqrt(diag(Vmatrix_pred))
y_pred_lwr <- yhat_test - diff
y_pred_upr <- yhat_test + diff
Dtest$yhat_test_lower <- y_pred_lwr
Dtest$yhat_test_upper <- y_pred_upr


## 3.4 Plotting model and forecast with prediction interval
ggplot(Dtrain, aes(x = year)) +
  geom_point(aes(y = total, color = "Observed")) +
  geom_line(aes(y = yhat, color = "Fitted"), size = 0.5) +
  geom_point(data = Dtest,
             aes(x = year, y = yhat, color = "Forecast"),
             size = 1) +
  geom_ribbon(data = Dtest,
              aes(x = year,
                  ymin = yhat_test_lower,
                  ymax = yhat_test_upper,
                  fill = "Prediction Interval"),
              alpha = 0.2) +
  scale_color_manual(
    name = "",
    values = c("Observed" = "black",
               "Fitted" = "black",
               "Forecast" = "red")
  ) +
  scale_fill_manual(
    name = "",
    values = c("Prediction Interval" = "red")
  ) +
  xlab("Time [monthly]") +
  ylab("Number of cars registered [millions]") +
  ggtitle("OLS model with forecast and prediction interval") +
  theme(plot.title = element_text(hjust = 0.5)) 


## 3.5 Comment on prediction ----
# Taking the more or less constant growth of the number of registered cars from 2022 to the end 
# of 2023 into account, the forecast seems a bit too optimistic. 


## 3.6 Inspecting residuals
Dtrain$residual <- e_OLS
ggplot(Dtrain, aes(x = year, y=residual)) +
  geom_point(color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size=0.8) +
  xlab("Time [monthly]") +
  ylab("Residual") +
  ggtitle("Residuals for the OLS model") +
  theme(plot.title = element_text(hjust = 0.5))

# The residuals seem to have mean approximately equal to zero 
# But the residuals does not seem to be independent from the time 

