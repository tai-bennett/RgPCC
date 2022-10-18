# RgPCR algorithms

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
  output <- list(gamma_hat, beta_hat, U)
  names(output) <- c("gamma_hat", "beta_hat", "U_pseudo")
  #print(V)
  #print(S)
  return(output)
}


### --------------------------------- Zou's RgPCR Algorithm
# --------- Zou's RgPCR with Lasso
RgPCR_lasso_adjusted_v2 <- function(X, Y, lambda, U_old) {
  #cat("\n", "dim S", dim(S))
  #cat("\n", "dim u", dim(svd(X)$u))
  p <- ncol(X)
  U <- svd(X)$u
  D <- svd(X)$d
  V <- svd(X)$v

  S <- reflection_matrix_v2(U_old, U)
  #cat("\n", "dim S", dim(S))
  #cat("\n", "dim u", dim(svd(X)$u))
  U <- U %*% S
  V <- V %*% S
  
  gamma_hat_ols <- t(Y) %*% U
  gamma_hat <- t(pos_part(abs(gamma_hat_ols) - (0.5 * lambda)*D^(-1)) * sign(gamma_hat_ols))
  s <- t(as.matrix(gamma_hat / D))
  #beta_hat <- V %*% (gamma_hat / D)
  output <- list(gamma_hat, U, D, V)
  names(output) <- c("gamma_hat", "U_pseudo", "D", "V_pseudo")
  #print(V)
  #print(S)
  return(output)
}


