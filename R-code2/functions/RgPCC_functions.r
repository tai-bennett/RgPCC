
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

# =================================================================
# RgPCC PREDICTION
# =================================================================
# to predict and calculate testing error of RgPCC model
# -----------------------------------------------------------------

RgPCC.predict <- function(X.tune, Y.tune, beta.hat) {
  p.hat <- matrix(0, nrow(X.tune), 1)
  for (j in seq_len(nrow(X.tune))) {
    p.hat[j, ] <- mylogistic(X.tune[j, ], as.matrix(beta.hat))
  }
  
  y.hat.tune <- I(p.hat > 0.5)
  MSE <- mean(I(y.hat.tune != Y.tune))
  SSE <- sum(I(y.hat.tune != Y.tune))
  output <- list(p.hat, MSE, SSE)
  names(output) <- c("p.hat", "MSE", "SSE")
  return(output)
}

get.MSECV <- function(fold.size, X.train, Y.train, lambda, tol_0, Sigma){
  folds <- sample(fold.size, nrow(X.train), replace = T)
  MSEparts <- matrix(0, 1, fold.size)
  
  for (j in c(1:fold.size)){ 
    # train for each fold
    RgPCC_lasso_results <- RgPCC_lasso_experimental_v2(
      X.train[which(folds != j),], Y.train[which(folds != j),], lambda, tol_0, Sigma
    )
    # use model to predict and find error
    #print(RgPCC.predict(X.train[which(folds != j),], Y.train[which(folds != j),], RgPCC_lasso_results$beta_hat)$MSE)
    MSEparts[1,j] <- RgPCC.predict(X.train[which(folds == j),], Y.train[which(folds == j),], RgPCC_lasso_results$beta_hat)$SSE
  }
  #print(MSEparts)
  output <- list(sum(MSEparts)/nrow(X.train))
  names(output) <- c("CVerror")
  return(output)
}

get.MSECValt <- function(fold.size, X.train, Y.train, lambda, tol_0, Sigma){
  folds <- sample(fold.size, nrow(X.train), replace = T)
  MSEparts <- matrix(0, 1, fold.size)
  
  for (k in c(1:fold.size)){ 
    # train for each fold
    RgPCC_lasso_results <- RgPCC_lasso_experimental_v3(
      X.train[which(folds != k),], Y.train[which(folds != k),], lambda, tol_0, Sigma
    )
    # use model to predict and find error
    #print(RgPCC.predict(X.train[which(folds != j),], Y.train[which(folds != j),], RgPCC_lasso_results$beta_hat)$MSE)
    MSEparts[1,k] <- RgPCC.predict(X.train[which(folds == k),], Y.train[which(folds == k),], RgPCC_lasso_results$beta_hat)$SSE
  }
  #print(MSEparts)
  output <- list(sum(MSEparts)/nrow(X.train))
  names(output) <- c("CVerror")
  return(output)
}




# =================================================================
# Collecting the results of the RgPCC
# =================================================================
# 
# -----------------------------------------------------------------

metrics <- function(RgPCC.results, N, X.test, Y.test, p.test) {
  
  # storage
  storage <- matrix(0, N, 7)
  for (i in 1:N){
    # get beta and p
    beta.hat <- RgPCC.results[[i]]$beta_hat
    p.hat <- RgPCC.predict(X.test, Y.test, beta.hat)$p.hat
    gamma.size <- RgPCC.results[[i]]$gamma_size
    # test error
    class.error <- RgPCC.predict(X.test, Y.test, beta.hat)$MSE
    # p info
    p.info <- get.my.metrics(p.hat, p.test)
	# time and iterations
	time <- RgPCC.results[[i]]$time
	iter <- RgPCC.results[[i]]$total_iter
    # store it
    storage[i,] <- c(p.info, class.error, gamma.size, time, iter)
  }
  # calc mean of each metric
  output <- colMeans(storage)
  names(output) <- c("one.norm", "two.norm", "EGKL", "class.error", "gamma.size", "time", "total_iter")
  return(output)
}

# =================================================================
# a way to get the tuning parameters
# =================================================================
# 
# -----------------------------------------------------------------


store.tuning.data <- function(sample.size, lambda.set, tune.data){
  return(cbind(sample.size, lambda.set, colMeans(tune.data)))
}

