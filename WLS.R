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
lambda = 0.9
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

## 4.4 ----
y_train <- Dtrain$total
X_train <- cbind(1, Dtrain$year)

### WLS ----
theta_WLS <- solve(t(X_train)%*%solve(SIGMA_local)%*%X_train)%*%(t(X_train)%*%solve(SIGMA_local)%*%y_train)
print(theta_WLS)

# Check with built in model:
WLS_lm <- lm(total ~ year, weights = weights, data = Dtrain)

summary(WLS_lm)

# Predict on train
yhat_train_WLS <- X_train%*%theta_WLS

# Make to dataframe

train_preds_WLS <- X_train %>% 
  cbind(yhat_train_WLS) %>% 
  as.data.frame() %>% 
  select(year=V2, y_hat=V3)%>% 
  mutate(set='train')


### OLS ----

theta_OLS <- solve(t(X_train)%*%solve(SIGMA_global)%*%X_train)%*%(t(X_train)%*%solve(SIGMA_global)%*%y_train)
print(theta_OLS)
yhat_train_OLS <- X_train %*% theta_OLS

OLS_lm <- lm(total ~ year, data = Dtrain)

summary(OLS_lm)

train_preds_OLS <- X_train %>% 
  cbind(yhat_train_OLS) %>% 
  as.data.frame() %>% 
  select(year=V2, y_hat=V3)%>% 
  mutate(set='train')

## 4.5: Test data ----

y_test <- Dtest$total
X_test <- cbind(1, Dtest$year)

# Make test predictions
yhat_test_WLS <- X_test%*%theta_WLS
yhat_test_OLS <- X_test%*%theta_OLS

preds_LM_WLS <- (predict(WLS_lm, newdata = Dtest, interval = 'confidence'))
lower <- preds_LM_WLS[1:nrow(Dtest),2]
upper <- preds_LM_WLS[1:nrow(Dtest),3]

test_preds_WLS <- X_test %>% 
  cbind(yhat_test_WLS) %>% 
  as.data.frame() %>% 
  select(year=V2, y_hat=V3)%>% 
  mutate(set='test') %>% 
  cbind(lower) %>% 
  cbind(upper)

preds_WLS <- D %>% 
  left_join(bind_rows(train_preds_WLS, test_preds_WLS), by = c('year'))

ggplot(data=preds_WLS, aes(x=time, y=total, colour = set)) +
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


# OLS 

preds_LM_OLS <- (predict(OLS_lm, newdata = Dtest, interval = 'confidence'))
lower_OLS <- preds_LM_OLS[1:nrow(Dtest),2]
upper_OLS <- preds_LM_OLS[1:nrow(Dtest),3]


test_preds_OLS <- X_test %>% 
  cbind(yhat_test_OLS) %>% 
  as.data.frame() %>% 
  select(year=V2, y_hat=V3) %>% 
  mutate(set='test') %>% 
  cbind(lower_OLS) %>% 
  cbind(upper_OLS)

preds_OLS <- D %>% 
  left_join(bind_rows(train_preds_OLS, test_preds_OLS), by = c('year'))


ggplot(data=preds_OLS, aes(x=time, y=total, colour = set)) +
  geom_point() +
  geom_line(aes(x=time, y=y_hat, col='Predictions'),) +
  geom_ribbon(aes(ymin = lower_OLS, ymax=upper_OLS), alpha=0.2) + 
  geom_vline(xintercept = teststart) +
  xlab("Time [Month]") +
  ylab("Predicted [line] vs. actual [points]") +
  ggtitle("Forecast for OLS") +
  scale_color_manual(name = "Info"
                     , values = c("Predictions" = "red",'test'='blue', 'train'='grey'))


preds_WLS <- preds_WLS %>% 
  mutate(residuals = total - y_hat)


ggplot(data=preds_WLS, aes(x=time, y=residuals)) + 
  geom_point()

