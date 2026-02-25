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

# 4. WLS ----

## 4.1 ----

var_cov_matrix <- function(n,weights=NA) {
  if (all(is.na(c(weights)))) {
    weights <- rep(1,n)
  }
  SIGMA <- diag(n)
  diag(SIGMA) <- 1/weights
  return (SIGMA)
}

n <- nrow(Dtrain)
lambda <- 0.9
weights <- lambda^((n-1):0)


SIGMA_local <- var_cov_matrix(n, weights)
# print lower right corner to check:
print(SIGMA_local[(n-5):n,(n-5):n])

SIGMA_global <- var_cov_matrix(n)
print(SIGMA_global[(n-5):n,(n-5):n])

# In local model, we weight points further away with a higher variance, thus giving less weight (because more uncertainty)
# In global model, all variances are 1, meaning no difference in certainty

## 4.2: Plot lambda weights ----

Dtrain <- Dtrain %>% 
  mutate(weights=weights)

ggplot(data= Dtrain, aes(x=time,y=weights)) +
  geom_bar(stat='identity') +
  xlab("Data point [idx]") +
  ylab(latex2exp::TeX("$\\lambda$")) +
  ggtitle("Weights vs. time, the most recent point has highest weight")


## 4.3: Sum of weights ----

# WLS (local trend)
Dtrain %>% 
  summarise(sum(weights))

# OLS (global trend)

n

# OLS has n as the sum of weights


y_train <- Dtrain$total
X_train <- cbind(1, Dtrain$year)
y_test <- Dtest$total
X_test <- cbind(1, Dtest$year)
preds_WLSs <- list()

n <- nrow(Dtrain)
lambdas <- c(0.61, 0.7, 0.8, 0.9, 0.99)

for (i in 1:length(lambdas)) {
  lambda <- lambdas[i]
  weights <- lambda^((n-1):0)
  
  SIGMA_local <- var_cov_matrix(n, weights)
  
  
  # Estimate theta using equation
  theta_WLS <- solve(t(X_train)%*%solve(SIGMA_local)%*%X_train)%*%(t(X_train)%*%solve(SIGMA_local)%*%y_train)
  
  # Check with built in model:
  WLS_lm <- lm(total ~ year, weights = weights, data = Dtrain)
  
  summary(WLS_lm)
  
  # Predict on train and test
  # Make a dataframe
  yhat_train_WLS <- X_train%*%theta_WLS
  yhat_test_WLS <- X_test%*%theta_WLS
  
  preds_LM_WLS <- (predict(WLS_lm, newdata = Dtest, interval = 'confidence'))
  lower <- preds_LM_WLS[1:nrow(Dtest),2]
  upper <- preds_LM_WLS[1:nrow(Dtest),3]
  
  
  train_preds_WLS <- X_train %>% 
    cbind(yhat_train_WLS) %>% 
    as.data.frame() %>% 
    select(year=V2, y_hat=V3)%>% 
    mutate(set='train',
           lambda = lambda)
  
  test_preds_WLS <- X_test %>% 
    cbind(yhat_test_WLS) %>% 
    as.data.frame() %>% 
    select(year=V2, y_hat=V3)%>% 
    mutate(set='test',
           lambda = lambda) %>% 
    cbind(lower) %>% 
    cbind(upper)
  
  
  preds_WLS <- D %>% 
    left_join(bind_rows(train_preds_WLS, test_preds_WLS), by = c('year'))
  
  preds_WLSs[[i]] <- preds_WLS
  
}


preds_WLS <- bind_rows(preds_WLSs)

ggplot(data=preds_WLS, aes(x=time, y=total, colour = set, group = lambda)) +
  facet_wrap(~lambda, ncol=1, labeller=label_both) + 
  geom_point() +
  geom_ribbon(aes(ymin = lower, ymax=upper), alpha=0.2) + 
  geom_line(aes(x=time, y=y_hat, col='Predictions'),) +
  geom_vline(xintercept = teststart) +
  xlab("Time [Month]") +
  ylab("Predicted [line] vs. actual [points]") +
  ylim(c(2.9,3.4)) +
  ggtitle("Forecast for WLS") +
  scale_color_manual(name = "Info"
                     , values = c("Predictions" = "red",'test'='blue', 'train'='grey'))
2
