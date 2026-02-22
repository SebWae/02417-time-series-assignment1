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


## 5.2 Update equations algorithm
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


## 5.3 Final RLS params
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


## 5.4 RLS with forgetting
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
lambda <- 0.7
thetas <- get_rls_params_w_forget(Dtrain$year, Dtrain$total, theta_0, R_0, lambda)

# Unpacking the parameter values
for (theta in thetas) {
  
}


