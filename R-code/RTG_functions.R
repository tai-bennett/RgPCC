# Functions necessary of RgPCC_lasso


## packages
library(MASS)
library(Rlab)
library(matlib)


## --------------------------------------------- various functions


#### logistic function
mylogistic <- function(x, beta) {
    exp(x %*% beta) / (1 + exp(x %*% beta))
    }



#### positive part function
pos_part <- function(x) {
  for (i in 1:length(x)) {
    if (x[i] < 0) {
        x[i] <- 0
        }
        }
  return(x)
    }



#### indicator function
I <- function(condition) {
    ifelse(condition, 1, 0)
    }



# get reflection
get_reflection <- function (v1, v2) {
  ifelse(
    sum(v1 * v2) > 0,
    1,
    -1
  )
}

# reflection matrix
reflection_matrix <- function (Sigma, data) {
  #print(dim(Sigma))
  v_true <- eigen(Sigma)$vector
  v_data <- svd(data)$v
  #printMatrix(v_data)
  p <- ncol(data)
  N <- nrow(data)
  s <- diag(rep(0, min(c(p,N))))
  for (j in 1:min(c(p,N))) {
    s[j,j] <- get_reflection(v_true[, j], v_data[, j])
  }
  return(s)
}



### --------------------------------- simulated data function
sim_data <- function(size, mean, Sigma, gamma, seed) {
# ---------------- set seed
set.seed(seed)
# ---------------- simulate predictor data
X <- as.matrix(mvrnorm(n = size, mean, Sigma))
# ---------------- simulate response data
# oriented svd decomposition
  S <- reflection_matrix(Sigma, X)
  U <- svd(X)$u %*% S
  V <- svd(X)$v %*% S
  D <- svd(X)$d
# assign responses
  mu <- U %*% gamma
  prob <- exp(mu)/(1+exp(mu))
  Y <- as.matrix(rbern(size, prob))
  beta <- t(V)%*%inv(diag(D))%*%as.matrix(gamma)
  output <-list(X, Y, prob, beta)
  names(output) <- c("X", "Y", "prob", "true_beta")
  return(output)
}



### --------------------------------- Zou's RgPCR Algorithm
# --------- Zou's RgPCR with Lasso
RgPCR_lasso <- function(X, Y, lambda) {
p <- ncol(X)
U <- svd(X)$u
D <- svd(X)$d
V <- svd(X)$v
gamma_hat_ols <- t(Y) %*% U
gamma_hat <- t(pos_part(abs(gamma_hat_ols) - (0.5 * lambda)*D^(-1)) * sign(gamma_hat_ols))
s <- t(as.matrix(gamma_hat / D))
beta_hat <- V %*% (gamma_hat / D)
output <- list(gamma_hat, beta_hat)
names(output) <- c("gamma_hat", "beta_hat")
return(output)
}



### --------------------------------- Zou's RgPCR Algorithm
# --------- Zou's RgPCR with Lasso
RgPCR_lasso_adjusted <- function(X, Y, lambda, Cov_true) {
  S <- reflection_matrix(Cov_true, X)
  #cat("\n", "dim S", dim(S))
  #cat("\n", "dim u", dim(svd(X)$u))
  p <- ncol(X)
  U <- svd(X)$u %*% S
  D <- svd(X)$d
  V <- svd(X)$v %*% S
  gamma_hat_ols <- t(Y) %*% U
  gamma_hat <- t(pos_part(abs(gamma_hat_ols) - (0.5 * lambda)*D^(-1)) * sign(gamma_hat_ols))
  s <- t(as.matrix(gamma_hat / D))
  beta_hat <- V %*% (gamma_hat / D)
  output <- list(gamma_hat, beta_hat)
  names(output) <- c("gamma_hat", "beta_hat")
  #print(V)
  #print(S)
  return(output)
}


get.my.metrics <- function(prob, true.prob){
  one.norm <- sum(abs(prob - true.prob))
  two.norm <- sum((prob - true.prob)^2)
  EGKL <- sum(true.prob * log(true.prob/prob))
  
  output <- c(one.norm, two.norm, EGKL)
  names(output) <- c("one.norm", "two.norm", "EGKL")
  return(output)
}

#training and testing (pca) logistic
logistic.process <- function(X.train, Y.train, p.train, X.test, Y.test, p.test){
  
  #training and testing logistic regression
  log_fit_data.frame <- as.data.frame(cbind(Y.train, X.train))
  log_fit_test_data.frame <- as.data.frame(cbind(Y.test, X.test))
  
  logistictime.start <- Sys.time()
  logistic_fit <- glm(V1 ~ . -1 , data=log_fit_data.frame, family=binomial(link="logit"))
  
  D <- svd(X.train)$d
  V <- svd(X.train)$v
  beta <- logistic_fit$coef
  
  gamma <- diag(D) %*% V %*% beta
  
  gamma.size <- sum(gamma != 0)
  #print(logistic_fit)
  
  logistictime.end <- Sys.time()
  #cat("Logistic time :", difftime(logistictime.end, logistictime.start, units = time_units), time_units, "\n")
  
  p.hat.train.log <- predict(logistic_fit, log_fit_data.frame[2:(p+1)], type = "response")
  p.hat.test.log <- predict(logistic_fit, log_fit_test_data.frame[2:(p+1)], type = "response")
  y.hat.train.log <- I(p.hat.train.log > 0.5)
  y.hat.test.log <- I(p.hat.test.log > 0.5)
  mymetrics.test.log <- get.my.metrics(p.hat.test.log, p.test)
  class.error.test.log <- mean(I(y.hat.test.log != Y.test))
  

  # PCA
  pca <- princomp(X.train, cor=F)
  cumsum.var <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
  num.pc <- min(which(cumsum.var > 0.90))
  #cat("\n", "num.pc = ", num.pc)
  pca.X <- pca$scores[,c(1:num.pc)]
  pca.YX <- as.data.frame(cbind(Y.train, pca.X))
  pca.logistic_fit <- glm(V1 ~ . -1 , data=pca.YX, family=binomial(link="logit"))
  
  gamma.size.pca <- sum(pca.logistic_fit$coeff != 0)
  pca.X.test <- as.data.frame(X.test %*% pca$loadings)
  pca.X.test <- pca.X.test[,c(1:num.pc)]
  
  p.hat.train.logpca <- predict(pca.logistic_fit, as.data.frame(pca.X), type = "response")
  p.hat.test.logpca <- predict(pca.logistic_fit, as.data.frame(pca.X.test), type = "response")
  y.hat.train.logpca <- I(p.hat.train.logpca > 0.5)
  y.hat.test.logpca <- I(p.hat.test.logpca > 0.5)
  mymetrics.test.logpca <- get.my.metrics(p.hat.test.logpca, p.test)
  class.error.test.logpca <- mean(I(y.hat.test.logpca != Y.test))
  
  
  
  
  
  
  
  output <- list(p.hat.test.log,
                 y.hat.test.log,
                 mymetrics.test.log,
                 class.error.test.log,
                 p.hat.test.logpca,
                 y.hat.test.logpca,
                 mymetrics.test.logpca,
                 class.error.test.logpca,
                 gamma.size,
                 gamma.size.pca)
  
  names(output) <- c("p.hat.test.log",
                     "y.hat.test.log",
                     "mymetrics.test.log",
                     "class.error.test.log",
                     "p.hat.test.logpca",
                     "y.hat.test.logpca",
                     "mymetrics.test.logpca",
                     "class.error.test.logpca",
                     "gamma.size.log",
                     "gamma.size.pcalog")
  
  return(output)
}

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





metrics <- function(RgPCC.results, N, X.test, Y.test, p.test) {
  
  # storage
  storage <- matrix(0, N, 5)
  for (i in 1:N){
    # get beta and p
    beta.hat <- RgPCC.results[[i]]$beta_hat
    p.hat <- RgPCC.predict(X.test, Y.test, beta.hat)$p.hat
    gamma.size <- RgPCC.results[[i]]$gamma_size
    # test error
    class.error <- RgPCC.predict(X.test, Y.test, beta.hat)$MSE
    # p info
    p.info <- get.my.metrics(p.hat, p.test)
    # store it
    storage[i,] <- c(p.info, class.error, gamma.size)
  }
  # calc mean of each metric
  output <- colMeans(storage)
  names(output) <- c("one.norm", "two.norm", "EGKL", "class.error", "gamma.size")
  return(output)
}

# ========================================================================
# GRAPHICS
# ========================================================================

store.tuning.data <- function(sample.size, lambda.set, tune.data){
  return(cbind(sample.size, lambda.set, colMeans(tune.data)))
}

tune.visuals <- function(AIC, BIC, MSE, pMSE, MSECV, save.name){
  
  pdf(paste("(", save.name,",", p, ")_graphtuningAIC.pdf", sep = ""))
  myplotAIC <- ggplot(data = AIC, aes(x = lambda, y = AIC, group = sample_size)) +
    geom_line(aes(color = factor(sample_size))) +
    geom_point(aes(color = factor(sample_size))) +
    labs(title = "Parameter Tuning AIC", subtitle = paste("(", save.name,",", p, ")"), color = "Sample Size", y = "AIC")
  
  print(myplotAIC)
  dev.off()
  
  pdf(paste("(", save.name,",", p, ")_graphtuningBIC.pdf", sep = ""))
  myplotBIC <- ggplot(data = BIC, aes(x = lambda, y = BIC, group = sample_size)) +
    geom_line(aes(color = factor(sample_size))) +
    geom_point(aes(color = factor(sample_size))) +
    labs(title = "Parameter Tuning BIC", subtitle = paste("(", save.name,",", p, ")"), color = "Sample Size", y = "BIC")
  
  print(myplotBIC)
  dev.off()
  
  pdf(paste("(", save.name,",", p, ")_graphtuningMSE.pdf", sep = ""))
  myplotMSE <- ggplot(data = MSE, aes(x = lambda, y = MSE, group = sample_size)) +
    geom_line(aes(color = factor(sample_size))) +
    geom_point(aes(color = factor(sample_size))) +
    labs(title = "Parameter Tuning MSE", subtitle = paste("(", save.name,",", p, ")"), color = "Sample Size", y = "MSE")
  
  print(myplotMSE)
  dev.off()
  
  pdf(paste("(", save.name,",", p, ")_graphtuningpMSE.pdf", sep = ""))
  myplotpMSE <- ggplot(data = pMSE, aes(x = lambda, y = pMSE, group = sample_size)) +
    geom_line(aes(color = factor(sample_size))) +
    geom_point(aes(color = factor(sample_size))) +
    labs(title = "Parameter Tuning pMSE", subtitle = paste("(", save.name,",", p, ")"), color = "Sample Size", y = "pMSE")
  
  print(myplotpMSE)
  dev.off()
  
  pdf(paste("(", save.name,",", p, ")_graphtuningMSECV.pdf", sep = ""))
  myplotMSECV <- ggplot(data = MSECV, aes(x = lambda, y = MSECV, group = sample_size)) +
    geom_line(aes(color = factor(sample_size))) +
    geom_point(aes(color = factor(sample_size))) +
    labs(title = "Parameter Tuning MSECV", subtitle = paste("(", save.name,",", p, ")"), color = "Sample Size", y = "MSECV")
  
  print(myplotMSECV)
  dev.off()
  
}

metric.visuals <- function(results, save.name){
  
  #results <- results[order(results$method),]
  #results[,-1] <- round(as.numeric(results[,-1]), 4)
  
  g <- tableGrob(results)
  h <- grid::convertHeight(sum(g$heights), "mm", TRUE)
  w <- grid::convertHeight(sum(g$widths), "mm", TRUE)
  ggplot2::ggsave(paste("(", save.name,",", p, ")_metrics.pdf", sep = ""), g, height = h, width = w+3, units = "mm")
}

