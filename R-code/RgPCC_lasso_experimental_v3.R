# The RgPCC Algorithm with Lasso

#region [imported packages]
## packages
library(MASS)
library(Rlab)
library(matlib)
#endregion

#region [RgPCC_lasso algo]
 # ----------------------------------- RgPCC with lasso
RgPCC_lasso_experimental_v3 <- function(X, Y, lambda, tolerance, print = TRUE, singular.adjust = 0.000001, loop_limit = 50) {
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

# ======================================================================
# LOOP START LOOP START LOOP START LOOP START LOOP START LOOP START 
# ======================================================================




  
# --------------------------------- loop but capped at 100 iterations
  print(error > tolerance)
  print(j < loop_limit)
while (error > tolerance & j < loop_limit) {
# --------------------------------- define prob vector
    # if first time, use X to calculate probs
    if (j == 1) {
        for (i in 1:nrow(X)) {
            prob[i, 1] <- mylogistic(X[i,], beta)
            if (is.nan(prob[i,1]) == TRUE){prob[i,1] <- 1}
                                        #cat("\n", "prob testing", prob[i,1])
            if (prob[i,1] > 1-singular.adjust){prob[i,1] <- 1-singular.adjust}
            if (prob[i,1] < singular.adjust){prob[i,1] <- singular.adjust}
        }
    } else {
    # usually we use U_tilde gamma to calculate probs
        U = inv(sqrt(W_old)) %*% U_old
        for (i in 1:nrow(U)) {
            prob[i, 1] <- mylogistic(U[i,], gamma)
            if (is.nan(prob[i,1]) == TRUE){prob[i,1] <- 1}
                                        #cat("\n", "prob testing", prob[i,1])
            if (prob[i,1] > 1-singular.adjust){prob[i,1] <- 1-singular.adjust}
            if (prob[i,1] < singular.adjust){prob[i,1] <- singular.adjust}
        }
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
  
  #printMatrix(prob_temp)
  if(nrow(X) < ncol(X)) {slim <- TRUE; error <- 0; if(print==TRUE){cat("\n", "thin SVD is TRUE (", nrow(X), ncol(X), ")" )}}
# --------------------------------- define pseudo data and response
  #if (nrow(X) - length(zero_list) >= ncol(X)) { # if not slim
  X_pseudo <- sqrt(W_temp) %*% X_temp
  #printMatrix(sqrt(W_temp))
    if (j == 1){
        z <- sqrt(W_temp) %*% X_temp %*% beta + inv(sqrt(W_temp)) %*% (Y_temp-prob_temp)
    } else {
        z <- sqrt(W_temp) %*% inv(sqrt(W_old)) %*% U_old %*% gamma + inv(sqrt(W_temp)) %*% (Y_temp-prob_temp)
    }
# --------------------------------- Call Zou's algo and define new gamma (beta)
    if (j == 1){
        coeff <- RgPCR_lasso_adjusted(X_pseudo, z, lambda, Sigma_temp)
        if (ncol(X_pseudo) == 2) {
            ggplot(data = data.frame(X_pseudo), aes(x = X1, y = X2)) + geom_point()
                                        #print(svd(X_pseudo)$v)
        }
    } else {
        coeff <- RgPCR_lasso_adjusted_v2(X_pseudo, z, lambda, U_old)
    }
# --------------------------------- test tolerance
  #beta_new <- coeff[[2]]
#  if (all(beta == matrix(0, nrow = ncol(X)))) {
#    error <- 99999
#    #cat("\n", "beta is zero")
#    }
#  else {
#    error <- sqrt(sum((beta - beta_new)^2)) / sqrt(sum(beta^2))
#    #cat("\n", "this is beta", beta)
#    #cat("\n", "this is gamma", coeff[[1]])
#    }
  gamma_new <- coeff[[1]]
    if (j > 1) {error = sqrt(sum((gamma - gamma_new)^2)) / sqrt(sum(gamma^2))}
# --------------------------------- update potential output
  output <- list(coeff[[1]], coeff[[2]], prob, error, j + 1 == loop_limit, j+1, slim, colSums(as.matrix(coeff[[1]] != 0)))
  #cat("\n", "this is gamma ", coeff[[1]], "\n", "this is beta", coeff[[2]])
# --------------------------------- update variables
  #beta <- beta_new
  gamma <- gamma_new
    if (j > 1){
        U_old <- coeff[[2]]
    } else {
        U_old <- coeff[[3]]
    }
    W_old <- W
  #cat("\n", "this is beta_new: ", beta_new)
  j <- j + 1
}


  
# ======================================================================
# LOOP END LOOP END LOOP END LOOP END LOOP END LOOP END LOOP END 
# ======================================================================

gamma_hat <- gamma
beta_hat <- coeff[[4]] %*% (gamma_hat / coeff[[3]])
output[[2]] <- beta_hat
  
  names(output) <-
  c("gamma_hat", "beta_hat", "prob_hat", "error", "max_interations?", "iterations", "slim?", "gamma_size")
  
  #if (length(zero_list) != 0) {cat("\n","max removed ", max(zero.lengths), " instances from X")}
  
  
  return(output)
}
#endregion
