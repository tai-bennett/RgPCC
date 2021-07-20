# The RgPCC Algorithm with Lasso

#region [imported packages]
## packages
library(MASS)
library(Rlab)
library(matlib)
#endregion

#region [RgPCC_lasso algo]
 # ----------------------------------- RgPCC with lasso
RgPCC_lasso_experimental_v2 <- function(X, Y, lambda, tolerance, print = TRUE, singular.adjust = 0.000001, loop_limit = 50) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Sigma <- t(X) %*% X
  error <- 99999
  #loop_limit <- 50
# ---------------------------------- initialize variables
  beta <- matrix(0, nrow = ncol(X))
  prob <- matrix(rep(0, nrow(X)), nrow = nrow(X))
  W <- matrix(rep(0, nrow(X) * nrow(X)), nrow = nrow(X))
  j <- 1
  zero.lengths <- c(0)
# --------------------------------- loop but capped at 100 iterations
while (error > tolerance & j < loop_limit) {
# --------------------------------- define prob vector
  for (i in 1:nrow(X)) {
    prob[i, 1] <- mylogistic(X[i,], beta)
    if (is.nan(prob[i,1]) == TRUE){prob[i,1] <- 1}
    #cat("\n", "prob testing", prob[i,1])
    if (prob[i,1] > 1-singular.adjust){prob[i,1] <- 1-singular.adjust}
    if (prob[i,1] < singular.adjust){prob[i,1] <- singular.adjust}
  }
# --------------------------------- define weight matrix
   for(i in 1:nrow(X)) {
     W[i, i] <- prob[i] * (1 - prob[i])
   }
  slim <- FALSE
  X_temp <- X
  Y_temp <- Y
  W_temp <- W
  prob_temp <- prob
  Sigma_temp <- Sigma
  zero_list <- which(diag(W) < 0.000001|is.nan(diag(W)))
  
  
  
  #if(length(zero_list) != 0 ){
  #W_temp <- as.matrix(W[-zero_list, -zero_list])
  
  #cat("\n", "diag of W_temp is ", diag(W_temp))
  #X_temp <- as.matrix(X[-zero_list,])
  #Y_temp <- as.matrix(Y[-zero_list])
  #prob_temp <- as.matrix(prob[-zero_list,])
  #Sigma_temp <- Sigma[-zero_list, -zero_list]
  #}
  
  #if(ncol(X_temp) == 1){X_temp <- t(X_temp)}
  #cat("\n", "dim of X_temp is", dim(X_temp))
  #cat("\n", "length of zero list", length(zero_list))
  
  #if (length(zero_list) != 0){zero.lengths <- c(zero.lengths, length(zero_list))}
  
  
  #fix NaN issue
  #prob_temp[is.nan(prob_temp)] <- 0.999999
  
  #printMatrix(prob_temp)
  if(nrow(X) < ncol(X)) {slim <- TRUE; error <- 0; if(print==TRUE){cat("\n", "thin SVD is TRUE (", nrow(X), ncol(X), ")" )}}
  #cat("\n",diag(W))
  #cat("\n",dim(X))
  #cat("\n",dim(Y))
  #cat("\n",dim(prob))
  
  #cat("\n", "W = ", diag(W))
# --------------------------------- define pseudo data and response
  #if (nrow(X) - length(zero_list) >= ncol(X)) { # if not slim
  X_pseudo <- sqrt(W_temp) %*% X_temp
  #printMatrix(sqrt(W_temp))
  z <- sqrt(W_temp) %*% X_temp %*% beta + inv(sqrt(W_temp)) %*% (Y_temp-prob_temp)
# --------------------------------- Call Zou's algo and define new beta
  coeff <- RgPCR_lasso_adjusted(X_pseudo, z, lambda, Sigma_temp)
  if (ncol(X_pseudo) == 2) {
    ggplot(data = data.frame(X_pseudo), aes(x = X1, y = X2)) + geom_point()
    #print(svd(X_pseudo)$v)
  }
  
# --------------------------------- test tolerance
  beta_new <- coeff[[2]]
  if (all(beta == matrix(0, nrow = ncol(X)))) {
    error <- 99999
    #cat("\n", "beta is zero")
    }
  else {
    error <- sqrt(sum((beta - beta_new)^2)) / sqrt(sum(beta^2))
    #cat("\n", "this is beta", beta)
    #cat("\n", "this is gamma", coeff[[1]])
    }
# --------------------------------- update potential output
  output <- list(coeff[[1]], coeff[[2]], prob, error, j + 1 == loop_limit, j+1, slim, colSums(as.matrix(coeff[[1]] != 0)))
  #cat("\n", "this is gamma ", coeff[[1]], "\n", "this is beta", coeff[[2]])
# --------------------------------- update variables
  beta <- beta_new
  #cat("\n", "this is beta_new: ", beta_new)
  j <- j + 1
  #if (j == loop_limit) {
  #  cat("\n", "Change in beta is ", error, " percent.")
  #  }
  #} end of "if not slim" statement
#cat("\n", "-")
}
  names(output) <-
  c("gamma_hat", "beta_hat", "prob_hat", "error", "max_interations?", "iterations", "slim?", "gamma_size")
  
  #if (length(zero_list) != 0) {cat("\n","max removed ", max(zero.lengths), " instances from X")}
  
  
  return(output)
}
#endregion