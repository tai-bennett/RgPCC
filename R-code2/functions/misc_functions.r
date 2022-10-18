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

# ========================================================================
# ========================================================================
# SOME SIMPLE FUNCTION
# ========================================================================
# ========================================================================



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

# reflection matrix
reflection_matrix_v2 <- function (U_old, data) {
  #print(dim(Sigma))
  v_true <- svd(U_old)$v
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

