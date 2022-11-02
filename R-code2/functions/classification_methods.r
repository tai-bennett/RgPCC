
## packages
library(MASS)
library(Rlab)
library(matlib)
library(xtable)
library(glmnet)
library(elasticnet)
library(caret)
library(grid)
library(gridExtra)
library(dplyr)
library(multcomp)
library(stargazer)
library(e1071)

source("misc_functions.r")

# =================================================================
# LOGISTIC REGRESSION and PCA LOGISTIC REGRESSION
# =================================================================
# 
# -----------------------------------------------------------------

#training and testing (pca) logistic
logistic.process <- function(X.train, Y.train, p.train, X.test, Y.test, p.test){
  
  #training and testing logistic regression
  log_fit_data.frame <- as.data.frame(cbind(Y.train, X.train))
  log_fit_test_data.frame <- as.data.frame(cbind(Y.test, X.test))
  
  logistictime.start <- Sys.time()
  logistic_fit <- glm(V1 ~ . -1 , data=log_fit_data.frame, family=binomial(link="logit"))
  logistictime.end <- Sys.time()
  time = logistictime.end - logistictime.start

  D <- svd(X.train)$d
  V <- svd(X.train)$v
  beta <- logistic_fit$coef
  
  gamma <- diag(D) %*% V %*% beta
  
  gamma.size <- sum(gamma != 0)
  #print(logistic_fit)
  
  #cat("Logistic time :", difftime(logistictime.end, logistictime.start, units = time_units), time_units, "\n")
  
  p.hat.train.log <- predict(logistic_fit, log_fit_data.frame[2:(p+1)], type = "response")
  p.hat.test.log <- predict(logistic_fit, log_fit_test_data.frame[2:(p+1)], type = "response")
  y.hat.train.log <- I(p.hat.train.log > 0.5)
  y.hat.test.log <- I(p.hat.test.log > 0.5)
  mymetrics.test.log <- get.my.metrics(p.hat.test.log, p.test)
  class.error.test.log <- mean(I(y.hat.test.log != Y.test))
  time.log <- time
  

  # PCA
  time.start <- Sys.time()
  pca <- princomp(X.train, cor=F)
  cumsum.var <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
  num.pc <- min(which(cumsum.var > 0.90))
  #cat("\n", "num.pc = ", num.pc)
  pca.X <- pca$scores[,c(1:num.pc)]    # 
  pca.YX <- as.data.frame(cbind(Y.train, pca.X))
  pca.logistic_fit <- glm(V1 ~ . -1 , data=pca.YX, family=binomial(link="logit"))
  time.stop <- Sys.time()
  time = time.stop - time.start
  gamma.size.pca <- sum(pca.logistic_fit$coeff != 0)
  pca.X.test <- as.data.frame(X.test %*% pca$loadings)
  pca.X.test <- pca.X.test[,c(1:num.pc)]
  
  p.hat.train.logpca <- predict(pca.logistic_fit, as.data.frame(pca.X), type = "response")
  p.hat.test.logpca <- predict(pca.logistic_fit, as.data.frame(pca.X.test), type = "response")
  y.hat.train.logpca <- I(p.hat.train.logpca > 0.5)
  y.hat.test.logpca <- I(p.hat.test.logpca > 0.5)
  mymetrics.test.logpca <- get.my.metrics(p.hat.test.logpca, p.test)
  class.error.test.logpca <- mean(I(y.hat.test.logpca != Y.test))
  time.logpca <- time
  
  
  
  
  
  
  
  output <- list(p.hat.test.log,
                 y.hat.test.log,
                 mymetrics.test.log,
                 class.error.test.log,
                 p.hat.test.logpca,
                 y.hat.test.logpca,
                 mymetrics.test.logpca,
                 class.error.test.logpca,
                 gamma.size,
                 gamma.size.pca,
				 time.log,
				 time.logpca)

  
  names(output) <- c("p.hat.test.log",
                     "y.hat.test.log",
                     "mymetrics.test.log",
                     "class.error.test.log",
                     "p.hat.test.logpca",
                     "y.hat.test.logpca",
                     "mymetrics.test.logpca",
                     "class.error.test.logpca",
                     "gamma.size.log",
                     "gamma.size.pcalog",
					 "time.log",
					 "time.logpca")

  return(output)
}

# =================================================================
# ELASTIC NET AND RIDGE REGRESSION
# =================================================================
# 
# -----------------------------------------------------------------

# training and testing for enet and ridge fit
# this trains enet with package 'caret'
enet.process.old <- function(X.train, Y.train, p.train, X.test, Y.test, p.test, lambda.set){
  cv_5 <- trainControl(method = "cv", number = 5)
  enet.grid <- expand.grid(alpha = seq(0,1,0.1),
                          lambda = lambda.set)
  ridge.grid <- expand.grid(alpha = 0,
                          lambda = lambda.set)
  
  enet.fit <- train(x = as.data.frame(X.train),
                    y = as.factor(Y.train),
                    method = "glmnet",
                    trControl = cv_5,
                    tuneGrid = enet.grid)

  p.hat.test.enet <- predict(enet.fit, newdata = as.data.frame(X.test), type = "prob")[[2]]
  y.hat.test.enet <- I(p.hat.test.enet > 0.5)
  mymetrics.test.enet <- get.my.metrics(p.hat.test.enet, p.test)
  class.error.test.enet <- mean(I(y.hat.test.enet != Y.test))
  gamma.size.enet <- sum(as.numeric(enet.fit$coeff) != 0)
  
  ridge.fit <- train(x = as.data.frame(X.train),
                    y = as.factor(Y.train),
                    method = "glmnet",
                    trControl = cv_5,
                    tuneGrid = ridge.grid)

  p.hat.test.ridge <- predict(ridge.fit, newdata = as.data.frame(X.test), type = "prob")[[2]]
  y.hat.test.ridge <- I(p.hat.test.ridge > 0.5)
  mymetrics.test.ridge <- get.my.metrics(p.hat.test.ridge, p.test)
  class.error.test.ridge <- mean(I(y.hat.test.ridge != Y.test))
  gamma.size.ridge <- sum(as.numeric(ridge.fit$coeff) != 0)                              
  
output <- list(p.hat.test.enet = p.hat.test.enet,
			   y.hat.test.enet = y.hat.test.enet,
               mymetrics.test.enet = mymetrics.test.enet,
               class.error.test.enet = class.error.test.enet,
			   gamma.size.enet = gamma.size.enet,
               p.hat.test.ridge = p.hat.test.ridge,
               y.hat.test.ridge = y.hat.test.ridge,
               mymetrics.test.ridge = mymetrics.test.ridge,
               class.error.test.ridge = class.error.test.ridge,
               gamma.size.ridge = gamma.size.ridge
               )
return(output)
}

enet.process <- function(X.train, Y.train, p.train, X.test, Y.test, p.test, lambda.set){
  time.start <- Sys.time()
  enet.fit <- cv.glmnet(X.train, Y.train, family = "binomial", type.measure = "class", intercept=FALSE)
  time.stop <- Sys.time()
  time.enet <- time.stop - time.start
  y.hat.test.enet <- predict(enet.fit, newx = X.test, type='class', s = enet.fit$lambda.min)
  p.hat.test.enet <- predict(enet.fit, newx = X.test, type='response', s = enet.fit$lambda.min)

  class.error.test.enet <- mean(I(as.integer(y.hat.test.enet) != Y.test))
  mymetrics.test.enet <- get.my.metrics(p.hat.test.enet, p.test)

  gamma = coef(enet.fit, s = enet.fit$lambda.min)
  gamma.size.enet <- sum(gamma != 0)

  time.start <- Sys.time()
  ridge.fit <- cv.glmnet(X.train, Y.train, family = "binomial", type.measure = "class", intercept=FALSE)
  time.stop <- Sys.time()
  time.ridge <- time.stop - time.start
  y.hat.test.ridge <- predict(ridge.fit, newx = X.test, type='class', s = ridge.fit$lambda.min)
  p.hat.test.ridge <- predict(ridge.fit, newx = X.test, type='response', s = ridge.fit$lambda.min)

  class.error.test.ridge <- mean(I(as.integer(y.hat.test.ridge) != Y.test))
  mymetrics.test.ridge <- get.my.metrics(p.hat.test.ridge, p.test)

  gamma = coef(enet.fit, s = enet.fit$lambda.min)
  gamma.size.ridge <- sum(gamma != 0)

output <- list(p.hat.test.enet = p.hat.test.enet,
			   y.hat.test.enet = y.hat.test.enet,
               mymetrics.test.enet = mymetrics.test.enet,
               class.error.test.enet = class.error.test.enet,
			   gamma.size.enet = gamma.size.enet,
               p.hat.test.ridge = p.hat.test.ridge,
               y.hat.test.ridge = y.hat.test.ridge,
               mymetrics.test.ridge = mymetrics.test.ridge,
               class.error.test.ridge = class.error.test.ridge,
               gamma.size.ridge = gamma.size.ridge,
			   time.enet = time.enet,
			   time.ridge = time.ridge
               )
return(output)
}
# =================================================================
# LASSO PENALIZED LOGISTIC REGRESSION
# =================================================================
# 
# -----------------------------------------------------------------

lasso.log.process <- function(
    X.train, 
    Y.train, 
    p.train, 
    X.test, 
    Y.test, 
    p.test, 
    lambda.set
    )
{

  time.start <- Sys.time()
  SVD.train <- svd(X.train)
  U.train <- SVD.train$u
  U.test <- X.test %*% SVD.train$v

  lasso.log.fit <- cv.glmnet(U.train, Y.train, family = "binomial", alpha = 1, type.measure = "class", intercept=FALSE)
  time.stop <- Sys.time()
  time.lasso.log <- time.stop - time.start
  y.hat.test <- predict(lasso.log.fit, newx = U.test, type='class', s = lasso.log.fit$lambda.min)
  p.hat.test <- predict(lasso.log.fit, newx = U.test, type='response', s = lasso.log.fit$lambda.min)

  class.error.test <- mean(I(as.integer(y.hat.test) != Y.test))
  mymetrics.test <- get.my.metrics(p.hat.test, p.test)

  gamma = coef(lasso.log.fit, s = lasso.log.fit$lambda.min)
  gamma.size <- sum(gamma != 0)

  output <- list(p.hat.test.lasso.log = p.hat.test,
		 y.hat.test.lasso.log = y.hat.test,
		 mymetrics.test.lasso.log = mymetrics.test,
		 class.error.lasso.log = class.error.test,
		 gamma.size.lasso.log = gamma.size,
		 time.lasso.log = time.lasso.log
		 )

  return(output)
}
# =================================================================
# GET.MY.LAMBDA
# =================================================================
# to get the optimal lambda for the RgPCC tuning process
# -----------------------------------------------------------------

get.my.lambda <- function(AIC, BIC, MSE, pMSE, MSECV, lambda.set){
  lambda.AIC <- lambda.set[which.min(colMeans(AIC))]
  lambda.BIC <- lambda.set[which.min(colMeans(BIC))]
  lambda.MSE <- lambda.set[which.min(colMeans(MSE))]
  lambda.pMSE <- lambda.set[which.min(colMeans(pMSE))]
  lambda.MSECV <- lambda.set[which.min(colMeans(MSECV))]
  
  lambda.AIC.index <- which.min(colMeans(AIC))
  lambda.BIC.index <- which.min(colMeans(BIC))
  lambda.MSE.index <- which.min(colMeans(MSE))
  lambda.pMSE.index <- which.min(colMeans(pMSE))
  lambda.MSECV.index <- which.min(colMeans(MSECV))
  
  output <- list(lambda.AIC, 
                 lambda.BIC, 
                 lambda.MSE, 
                 lambda.pMSE, 
                 lambda.MSECV, 
                 lambda.AIC.index, 
                 lambda.BIC.index, 
                 lambda.MSE.index, 
                 lambda.pMSE.index, 
                 lambda.MSECV.index)
  
  names(output) <- c("AIC", "BIC", "MSE", "pMSE", "MSECV", "AIC.index", "BIC.index", "MSE.index", "pMSE.index", "MSECV.index")
  return(output)
}



get.my.metrics <- function(prob, true.prob){
  N <- length(prob)
  one.norm <- sum(abs(prob - true.prob))/N
  two.norm <- sum((prob - true.prob)^2)/N
  EGKL <- sum(true.prob * log(true.prob/prob))/N
  
  output <- c(one.norm, two.norm, EGKL)
  names(output) <- c("one.norm", "two.norm", "EGKL")
  return(output)
}

