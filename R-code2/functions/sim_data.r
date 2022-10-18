
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


# simulate data function


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

