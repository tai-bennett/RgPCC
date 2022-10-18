# title: script1.r
# author: Duncan Bennett
# description: Compare the performance of RgPCC and penalized logistic
# regression via glmnet. The metrics used should be the same as the previous
# experiments. This calls experiment1.r which gives us the function that
# controls the experiment (similar to RgPCC_lasso_simulated_experiment.r)

source('RTG_functions.R')
source('experiment1.r')

# modifiable parameters
# ------------------------------------------------------------------------------
n <- 100 # sample size
p <- 15 # dimension of data
mu <- rep(0, p)# mean
rho <- 0.8 # covariance factor
gamma <- c(0, 0, 0, 20, 10, 5, rep(0, p-6))
seed <- 1234
# ------------------------------------------------------------------------------

# unmodifiable parameters
# ------------------------------------------------------------------------------
Sigma <- matrix(c(1:(p*p)), nrow = p); 
for(i in 1:p) {
  for (j in 1:p) {
    Sigma[i,j] <- (rho)^abs(i - j)
  }
}
# ------------------------------------------------------------------------------

# script
# make data
data.train <- sim_data(
  n,
  rep(0, p), 
  Sigma,
  gamma,
  seed
)

data.test <- sim_data(
  5*n,
  rep(0, p), 
  Sigma,
  gamma,
  seed+100
)
SVD.train = svd(data.train$X)
U.train = SVD.train$u
U.test = data.test$X %*% SVD.train$v

model <- lasso.log.process(
  U.train, 
  data.train$Y, 
  data.train$prob, 
  U.test,
  data.test$Y, 
  data.test$prob
)

model$lambda.min

# predict and calculate metrics
data.test.yhat <- predict(model, newx = U.test, type='class', s = model$lambda.min)
data.test.phat <- predict(model, newx = U.test, type='response', s = model$lambda.min)

misclassification_error = mean(I(as.integer(data.test.yhat) != data.test$Y))
prob_errors = get.my.metrics(data.test$prob, data.test.phat)

gamma = coef(model, s = model$lambda.min)
