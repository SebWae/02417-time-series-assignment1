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


# 5. Recursive estimation and optimization of lambda ----
## 5.1 Update equations by hand ----
# See written notes!
my_vect = c(1,1)
my_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
print(my_matrix)
mat_inv = solve(my_matrix)
print(mat_inv)
print(my_vect %*% t(my_vect))
print(2 - matrix(c(1), nrow = 1, ncol = 1))


## 5.2 Update equations algorithm ----
get_rls_params <- function(x_vals, y_vals, theta_0, R_0) {
  # Number of data points
  N <- length(x_vals)  
  
  # Lists to contain thetas and R-matrices  
  thetas <- vector("list", N)
  Rs <- vector("list", N)
  
  # Initializing theta_0 and R_0
  theta <- theta_0
  R <- R_0
  
  for (i in 1: N) {
    # Retrieving previous values
    theta_prev <- theta
    R_prev <- R
    
    # Retrieving data
    x_val <- x_vals[[i]]
    y_val <- y_vals[[i]]
    
    # Constructing design matrix
    X <- matrix(c(1, x_val), ncol = 1)
    
    # Computing next R matrix
    R <- R_prev + X %*% t(X)
    
    # Computing next parameter values
    #print(solve(R) %*% X)
    #print((y_val - t(X) %*% theta_prev))
    #print(solve(R) %*% X * (y_val - t(X) %*% theta_prev))
    theta <- theta_prev + (solve(R) %*% X) * as.numeric(y_val - t(X) %*% theta_prev)
    
    # Appending new values to lists
    thetas[[i]] <- theta
    Rs[[i]] <- R
  }
  return (thetas)
}

# Calling get_rls_params function
theta_0 <- c(0,0)
R_0 <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
thetas <- get_rls_params(Dtrain$year[1:3], Dtrain$total[1:3], theta_0, R_0)
print(thetas)

# The term Y_t - X_t^T @ theta_{t-1} is the error using the previous parameters
# If this term is large, the parameters will be updated by a larger magnitude
# If the error is 0, the parameters will not be updated at all 


## 5.3 Final RLS params ----
thetas <- get_rls_params(Dtrain$year, Dtrain$total, theta_0, R_0)
print(last(thetas))

# For t=N the parameters are (-0.05831943, 0.00156796)
# The OLS parameters were (-110.35542810, 0.05614456 
# They are not really close, the OLS intercept is much lower than the RLS intercept
# and thus the slope is also more steep for OLS

# Trying to modify the initial values such that the parameter estimates become more similar
theta_0 <- c(1,1)
R_0 <- matrix(c(0.000001, 0, 0, 0.000001), nrow = 2, ncol = 2)
thetas <- get_rls_params(Dtrain$year, Dtrain$total, theta_0, R_0)
print(last(thetas))

# When changing R_0 to a diagonal matrix with the diagonal entries equal to 1*10^{-6} 
# the parameter estimates become equal to (-108.30708807, 0.05513101)
# which is similar to the OLS estimates
# changing theta_0 to c(1,1) does not really make a difference
# choosing theta_0 = c(1,1) would be based on an initial belief in a positive correlation
# between time and the total number of registered vehicles 
# and an assumption that there are not 0 registered vehicles at time t=0
# changing R has a larger impact 


## 5.4 RLS with forgetting ----
get_rls_params_w_forget <- function(x_vals, y_vals, theta_0, R_0, lambda) {
  # Number of data points
  N <- length(x_vals)  
  
  # Lists to contain thetas and R-matrices  
  thetas <- vector("list", N)
  Rs <- vector("list", N)
  
  # Initializing theta_0 and R_0
  theta <- theta_0
  R <- R_0
  
  for (i in 1: N) {
    # Retrieving previous values
    theta_prev <- theta
    R_prev <- R
    
    # Retrieving data
    x_val <- x_vals[[i]]
    y_val <- y_vals[[i]]
    
    # Constructing design matrix
    X <- matrix(c(1, x_val), ncol = 1)
    
    # Computing next R matrix
    R <- lambda * R_prev + X %*% t(X)
    
    # Computing next parameter values
    #print(solve(R) %*% X)
    #print((y_val - t(X) %*% theta_prev))
    #print(solve(R) %*% X * (y_val - t(X) %*% theta_prev))
    theta <- theta_prev + (solve(R) %*% X) * as.numeric(y_val - t(X) %*% theta_prev)
    
    # Appending new values to lists
    thetas[[i]] <- theta
    Rs[[i]] <- R
  }
  return (thetas)
}

theta_0 <- c(0,0)
R_0 <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
thetas <- get_rls_params_w_forget(Dtrain$year, Dtrain$total, theta_0, R_0, 0.7)
print(thetas)
# Unpacking the parameter values
theta_0s <- c()
theta_1s <- c()

for (theta in thetas) {
  # Retrieving theta values
  theta_0 <- theta[1]
  theta_1 <- theta[2]
  
  # Appending to lists
  theta_0s <- append(theta_0s, theta_0)
  theta_1s <- append(theta_1s, theta_1)
}

Dtrain$theta_0_l07 <- theta_0s
Dtrain$theta_1_l07 <- theta_1s

# Repeating with lambda = 0.99
theta_0 <- c(0,0)
R_0 <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
thetas <- get_rls_params_w_forget(Dtrain$year, Dtrain$total, theta_0, R_0, 0.99)
print(thetas)
# Unpacking the parameter values
theta_0s <- c()
theta_1s <- c()

for (theta in thetas) {
  # Retrieving theta values
  theta_0 <- theta[1]
  theta_1 <- theta[2]
  
  # Appending to lists
  theta_0s <- append(theta_0s, theta_0)
  theta_1s <- append(theta_1s, theta_1)
}

Dtrain$theta_0_l099 <- theta_0s
Dtrain$theta_1_l099 <- theta_1s

# Plot for each parameter
ggplot(data=Dtrain[-(1:20), ], aes(x=time)) +
  geom_line(aes(y=theta_0_l07, col="lambda=0.7")) +
  geom_line(aes(y=theta_0_l099, col="lambda=0.99")) +
  xlab("Time [monthly]") +
  ylab(expression(theta[0])) +
  ggtitle(expression(paste("Values of ", theta[0]))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(
    name = "",
    values = c("lambda=0.7" = "red",
               "lambda=0.99" = "blue")
  ) 

ggplot(data=Dtrain[-(1:20), ], aes(x=time)) +
  geom_line(aes(y=theta_1_l07, col="lambda=0.7")) +
  geom_line(aes(y=theta_1_l099, col="lambda=0.99")) +
  xlab("Time [monthly]") +
  ylab(expression(theta[1])) +
  ggtitle(expression(paste("Values of ", theta[1]))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(
    name = "",
    values = c("lambda=0.7" = "red",
               "lambda=0.99" = "blue")
  ) 


# The parameter values seem to be stable at around 0 for lambda = 0.99
# More fluctuating parameter values for lambda = 0.7, especially for theta_0
# The parameter values does not seem to stabilize at some point for lambda = 0.7
# But they take some time to start fluctuating 


## 5.5 Predictions and residuals ----
N <- length(Dtrain$theta_0_l07)

# One-step predictions for lambda = 0.7
yhat_l07 <- c()
for (i in 1: N) {
  x_val <- Dtrain$year[[i]]
  if (i == 1) {
    theta_0 <- 0
    theta_1 <- 0
  }
  else {
    theta_0 <- Dtrain$theta_0_l07[[i-1]]
    theta_1 <- Dtrain$theta_1_l07[[i-1]]
  }
  pred <- theta_0 + x_val * theta_1
  yhat_l07 <- append(yhat_l07, pred)
}

# Adding predictions to Dtrain dataframe
Dtrain$yhat_l07 <- yhat_l07


# One-step predictions for lambda = 0.99
yhat_l099 <- c()
for (i in 1: N) {
  x_val <- Dtrain$year[[i]]
  if (i == 1) {
    theta_0 <- 0
    theta_1 <- 0
  }
  else {
    theta_0 <- Dtrain$theta_0_l099[[i-1]]
    theta_1 <- Dtrain$theta_1_l099[[i-1]]
  }
  pred <- theta_0 + x_val * theta_1
  yhat_l099 <- append(yhat_l099, pred)
}

# Adding predictions to Dtrain dataframe
Dtrain$yhat_l099 <- yhat_l099


# Computing residuals
Dtrain$res_l07 <- Dtrain$yhat_l07 - Dtrain$total
Dtrain$res_l099 <- Dtrain$yhat_l099 - Dtrain$total


# Plotting residuals
ggplot(data=Dtrain[-(1:4), ], aes(x=time)) +
  geom_line(aes(y=res_l07, col="lambda=0.7")) +
  geom_line(aes(y=res_l099, col="lambda=0.99")) +
  xlab("Time [monthly]") +
  ylab("Residual") +
  ggtitle("Residuals for one-step predictions") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(
    name = "",
    values = c("lambda=0.7" = "red",
               "lambda=0.99" = "blue")
  ) 

# The residuals do not have mean 0
# especially for lambda=0.99, where the residuals are negative at all times.
# For lambda=0.7, the residuals seem to approach 0, which indicates that the predictions become 
# more accurate over time. 
# As for the OLS residuals, the RLS residuals are not independent from the time
# The residuals look like a non-random time series.


## 5.6 Horizons and RMSE ----
N <- length(Dtrain$year)
horizons <- seq(from=1, to=12, by=1)   
lambdas <- seq(from=0.5, to=0.99, by=0.01)

# Matrix to store RMSEs
RMSE_matrix <- matrix(NA, nrow = length(lambdas), ncol = length(horizons))

for (h in horizons) {
  RMSEs <- c()
  for (lambda in lambdas) {
    # Get parameter estimates based on lambda
    theta_0 <- c(0,0)
    R_0 <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
    thetas <- get_rls_params_w_forget(Dtrain$year, Dtrain$total, theta_0, R_0, lambda)
    
    # Unpacking the parameter values
    theta_0s <- c()
    theta_1s <- c()
    
    for (theta in thetas) {
      # Retrieving theta values
      theta_0 <- theta[1]
      theta_1 <- theta[2]
      
      # Appending to lists
      theta_0s <- append(theta_0s, theta_0)
      theta_1s <- append(theta_1s, theta_1)
    }
    
    # Making predictions
    preds <- c()
    for (i in 1:N) {
      x_val <- Dtrain$year[[i]]
      if (i > h) {
        theta_0 <- theta_0s[[i-h]]
        theta_1 <- theta_1s[[i-h]]
        pred <- theta_0 + x_val * theta_1
        preds <- append(preds, pred)
      }
    }
    
    # Computing RMSE
    true_vals <- Dtrain$total[(h+1):N]
    residuals <- (preds - true_vals)**2
    RMSE = (1 / (N-h)) * sum(residuals)
    RMSEs <- append(RMSEs, RMSE)
  }
  RMSE_matrix[, h] <- RMSEs
}

# Plotting RMSEs against lambda values
# Set margins and remove grey background
df <- data.frame(lambdas, RMSE_matrix)

df_long <- df %>%
  pivot_longer(
    cols = -lambdas,
    names_to = "horizon",
    values_to = "RMSE"
  )

# Rename factor levels
df_long$horizon <- factor(df_long$horizon, levels = paste0("X", 1:12), labels = paste0("k=", 1:12))

# Make sure the color vector names match the new factor levels
colors <- c(
  "k=1"="#08306b", "k=2"="#08519c", "k=3"="#2171b5", 
  "k=4"="#4292c6", "k=5"="#6baed6", "k=6"="#9ecae1",
  "k=7"="#c6dbef", "k=8"="#dadaeb", "k=9"="#b2e2e2",
  "k=10"="#66c2a4", "k=11"="#2ca25f", "k=12"="#006d2c"
)

ggplot(df_long, aes(x = lambdas, y = RMSE, color = horizon)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  xlab(expression(lambda)) +
  ggtitle("RMSE for different horizons") +
  theme(plot.title = element_text(hjust = 0.5))

# The RMSE is proportional to value of lambda and the horizon k.
# Usually lambda is within the range (0.95, 0.999)
# A lower lambda value implies quicker forgetting 
# lambda = 1 corresponds to a static system in which the system parameters do not change over time (equal weight to all past data)
# It might be a good idea to look for an elbow point in the plot
# i.e., a point where the RMSE is not too high and the lambda is close to its typical range
# aka. a trade-off between low RMSE and a high lambda value
# lambda = 0.9 would be reasonable for short-term horizons
# lambda = 0.8 would be more appropriate for longer horizons (k=8 to k=12)
# So the choice of lambda depends on the horizon.


## 5.7 Test predictions ----
N_test <- length(Dtest$year)
lambdas <- seq(from=0.8, to=0.9, by=0.1/N_test)
preds <- c()

for (i in 1:N_test) {
  # Estimate parameters with lambda value
  lambda <- lambdas[[i]]
  theta_0 <- c(0,0)
  R_0 <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
  thetas <- get_rls_params_w_forget(Dtrain$year, Dtrain$total, theta_0, R_0, lambda)
  
  # Unpacking estimates
  theta_N <- thetas[[length(thetas)]]
  theta_0 <- theta_N[1, 1]
  theta_1 <- theta_N[2, 1]
  
  # Making prediction
  x_val <- Dtest$year[[i]]
  pred <- theta_0 + theta_1 * x_val
  preds <- append(preds, pred)
}

# Add predictions to test dataframe
Dtest$yhat_RLS <- preds


# Plotting training data and predictions
ggplot(Dtrain, aes(x = year)) +
  geom_point(aes(y = total, color = "Observed")) +
  geom_point(data = Dtest,
             aes(x = year, y = yhat_RLS, color = "Forecast")) +
  scale_color_manual(
    name = "",
    values = c("Observed" = "black",
               "Forecast" = "red")
  ) +
  xlab("Time [monthly]") +
  ylab("Number of cars registered [millions]") +
  ggtitle("RLS training data and predictions") +
  theme(plot.title = element_text(hjust = 0.5)) 


## 5.8 Reflections on time adaptive models ----
### Overfitting vs. underfitting ----
# For short horizons and low lambda values one might risk overfitting to the most recent data
# For longer horizons and lambda values close to 1, the model is more likely to underfit to new data
# and still be highly affected by the past.


### Creating time dependent test sets ----
# We cannot split the data in an arbitrary way, e.g., randomly sample 20% of the dataset.
# It must be sorted according to the time stamps. 
# If the data is recorded periodically, e.g., one record per day, 
# and the data follows some seasonal pattern 
# one might be careful when choosing the cutoff between the train and test set. 
# In most scenarios, it would be desirable to record data points with the same frequency across the 
# training and test set. Otherwise it will become more tricky to generalize to the unknown test set. 


### The role of recursive estimation in creating time dependent test sets ----
# If the training data follows some seasonal or periodic pattern, it might be straightforward 
# for RLS to take this into account, e.g., by letting the horizon be equal to 7 for a weekly pattern
# such that the parameter estimates from the same weekday last week are used to make predictions in this week


### Other techniques for time adaptive estimation ----
# Maybe some nested time series model can be applied in which the parameter values are being forecasted 
# by one or multiple separate model(s). 


### Additional thoughts ----
# Time adaptive models might be more susceptible to adversarial attacks since new malicious data
# will affect the model more compared to models that are not time adaptive.


