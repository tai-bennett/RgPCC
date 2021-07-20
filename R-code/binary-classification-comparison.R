# 5 fold cross validation and comparision of binary classificaiton methods



bin.method.comp <- function(X, Y, lambda.set, fold.size = 5, myseed, name, wd = getwd(), track = FALSE, print = FALSE) {
  original.wd <- getwd()
  library(glmnet)
  library(elasticnet)
  library(caret)
  library(grid)
  library(gridExtra)
  library(dplyr)
  library(multcomp)
  
  #source("D:/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/RgPCC_lasso_experimental.R")
  #source("D:/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/RTG_functions.R")
  
  source('RgPCC_lasso_experimental_v2.R')
  source('RTG_functions.R')
  setwd(wd)
  
  #source("C:/Users/benne/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/RgPCC_lasso_experimental.R")
  #source("C:/Users/benne/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/RTG_functions.R")
  
  time_units <- "sec"
  
  colnames(Y) <- "class"
  
  X <- scale(X, center = TRUE, scale = FALSE) # center the data
  
  YX <- cbind(Y,X)
  
  Y <- as.matrix(Y)
  
  # taken out, need to do PCA after folds
  #pca <- princomp(X[-which(folds == j)], cor=F)
  
  #cumsum.var <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
  #num.pc <- min(which(cumsum.var > 0.90))
  #cat("\n", "num.pc = ", num.pc)
  #pca.comp <- pca$scores[,c(1:num.pc)]
  #pca.YX <- as.data.frame(cbind(Y, pca.comp))

  #print(summary(pca))
# ----------------------------------------------------------------------------------- initializing for 5-fold cv loop
  
  set.seed(myseed)
  folds <- sample(fold.size, nrow(X), replace = T)
  
  N <- nrow(X)
  p <- ncol(X)
  
  beta.matrix <- matrix(0, p, fold.size)
  
  # empty matrix to hold results
  result.size <- length(lambda.set)
  results.LOO <- data.frame(
    lambda = rep(0, result.size),
    LOO_AIC = rep(0,result.size),
    LOO_BIC = rep(0,result.size),
    stringsAsFactors = FALSE)
  
  results.all <- data.frame(
    lambda = rep(0, result.size),
    AIC = rep(0,result.size),
    BIC = rep(0,result.size),
    stringsAsFactors = FALSE)
  
  results.AIC.LOO <- matrix(0, length(lambda.set), fold.size)
  results.BIC.LOO <- matrix(0, length(lambda.set), fold.size)
  
  # storing errors for comparison
  error.results <- data.frame(
    cverror = rep(0, 5*fold.size),
    method = rep(0, 5*fold.size)
  )  
  loglike.results <- data.frame(
    cverror = rep(0, 5*fold.size),
    method = rep(0, 5*fold.size)
  )


# ----------------------------------------------------------------------------------- cv loop

cv_5 <- trainControl(method = "cv", number = 5)
ridge.grid <- expand.grid(alpha = 0, 
                          lambda = lambda.set)


lambda.hat.AIC <- rep(0, fold.size) # indices for the optimal lambda
lambda.hat.BIC <- rep(0, fold.size)

y.hat.AIC <- rep(0, N) # vector of predicted class based on LOO models
y.hat.BIC <- rep(0, N)
y.hat.log <- rep(0, N)
y.hat.pcalog <- rep(0, N)
y.hat.enet <- rep(0, N)

AIC.data <- data.frame(row.names=1:4)
BIC.data <- data.frame(row.names=1:4)

for (j in c(1:fold.size)){ # cycle through each LOO
  
  AIC.vec <- rep(0, length(lambda.set))# AIC values for data-jth predictor over all lambdas
  BIC.vec <- rep(0, length(lambda.set))
  beta.matrix <- matrix(0, p, length(lambda.set))
  
  # RgPCC_lasso ------------------------------------------------------------------------------------------------------------
  
  for (i in seq(1, length(lambda.set))) { #through each lambda for TUNING

    if (i == 1) {start_time <- Sys.time(); 
    if (i == 1 & track==TRUE)cat("\n", "fold", j, "\n", "|")}
    if (0 == i%%1 & track==TRUE) {cat("=")}
    if (0 == i%%10 & track==TRUE) {cat("|")}
    modelRgPCC.LOO <- RgPCC_lasso_experimental_v2(X[-which(folds == j),], Y[-which(folds == j),], lambda.set[i], 0.1, print)
    beta.matrix[, i] <- modelRgPCC.LOO$beta_hat
    p.hat <- modelRgPCC.LOO$prob_hat
    loglikelihood <- sum(Y[-which(folds == j),] * log(p.hat) + (1-Y[-which(folds == j),]) * log(1-p.hat))
    gamma_size <- colSums(as.matrix(modelRgPCC.LOO$gamma_hat) != 0)
    
    AIC.vec[i] <- 2 * gamma_size - 2 * loglikelihood
    BIC.vec[i] <- gamma_size * log(N) - 2 * loglikelihood
    
    #results.LOO$lambda[i] <- lambda.set[i]
    #results.LOO$AIC[i] <- 2 * gamma_size - 2 * loglikelihood
    #results.LOO$BIC[i] <- gamma_size * log(N) - 2 * loglikelihood
    
  }
  
  AIC.data <- rbind(AIC.data, cbind(AIC.vec, lambda.set, j))
  BIC.data <- rbind(BIC.data, cbind(AIC.vec, lambda.set, j))
  
  #if(print==TRUE){cat("\n", "AIC", AIC.vec)}
  #if(print==TRUE){cat("\n", "BIC", BIC.vec)}
  lambda.hat.AIC[j] <- which.min(AIC.vec)
  lambda.hat.BIC[j] <- which.min(BIC.vec)
  
  p.hat.AIC <- mylogistic(X[which(folds == j),],beta.matrix[, lambda.hat.AIC[j]])
  p.hat.BIC <- mylogistic(X[which(folds == j),],beta.matrix[, lambda.hat.BIC[j]])
  loglike.results[j,] <- c(sum(Y[which(folds == j),] * log(p.hat.AIC) + (1-Y[which(folds == j),]) * log(1-p.hat.AIC)),1)
  loglike.results[j+5,] <- c(sum(Y[which(folds == j),] * log(p.hat.BIC) + (1-Y[which(folds == j),]) * log(1-p.hat.BIC)),2)
  
  
  y.hat.AIC[which(folds == j)] <- I(mylogistic(X[which(folds == j),],beta.matrix[, lambda.hat.AIC[j]]) > 0.5)
  y.hat.BIC[which(folds == j)] <- I(mylogistic(X[which(folds == j),],beta.matrix[, lambda.hat.BIC[j]]) > 0.5)
  
  error.results[j,] <- c(
    mean(I(y.hat.AIC[which(folds == j)] != t(Y)[which(folds == j)])),
    1)
  
  error.results[j+5,] <- c(
    mean(I(y.hat.BIC[which(folds == j)] != t(Y)[which(folds == j)])),
    2)
  
  # logistic regression ----------------------------------------------------------------------------------------------------
  
  logistic_fit <- glm(class ~ . -1 , data=as.data.frame(YX[-which(folds == j),]), family=binomial(link="logit"))
  
  p.hat.log <- as.numeric(predict(logistic_fit, as.data.frame(YX[which(folds == j),]), type = "response"))
  loglike.results[j+10,] <- c(sum(Y[which(folds == j),] * log(p.hat.log) + (1-Y[which(folds == j),]) * log(1-p.hat.log)),3)
  
  y.hat.log[which(folds == j)] <- I(as.numeric(predict(logistic_fit, as.data.frame(YX[which(folds == j),]), type = "response")) > 0.5)
  error.results[j+10,] <- c(
    mean(I(y.hat.log[which(folds == j)] != t(Y)[which(folds == j)])),
    3)
  
  # PCA logistic regression ------------------------------------------------------------------------------------------------
  
  # pca data frame
  pca <- princomp(X[-which(folds == j),], cor=F)
  
  cumsum.var <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
  num.pc <- min(which(cumsum.var > 0.90))
  cat("\n", "num.pc = ", num.pc)
  pca.comp <- pca$scores[,c(1:num.pc)]
  pca.YX <- as.data.frame(cbind(Y[-which(folds == j),], pca.comp))
  
  pca.logistic_fit <- glm(V1 ~ . -1 , data=pca.YX, family=binomial(link="logit"))
  
  #cat("\n", "Xfold dim is ", dim(X[which(folds == j),]))
  #cat("\n", "V dim is ", dim(pca$loadings))
  
  jfolddata <- as.data.frame(X[which(folds == j),] %*% pca$loadings)
  
  # convert testing data into eigenbasis and truncate by num.pc
  if(num.pc != 1){
    p.hat.pcalog <- as.numeric(predict(pca.logistic_fit, jfolddata[,c(1:num.pc)], type = "response"))
    loglike.results[j+15,] <- c(sum(Y[which(folds == j),] * log(p.hat.pcalog) + (1-Y[which(folds == j),]) * log(1-p.hat.pcalog)),4)
    y.hat.pcalog[which(folds == j)] <- I(as.numeric(predict(pca.logistic_fit, jfolddata[,c(1:num.pc)], type = "response")) > 0.5)
  }
  
  if(num.pc == 1){
    p.hat.pcalog <- as.numeric(predict(pca.logistic_fit, as.data.frame(t(jfolddata[,1])), type = "response"))
    loglike.results[j+15,] <- c(sum(Y[which(folds == j),] * log(p.hat.pcalog) + (1-Y[which(folds == j),]) * log(1-p.hat.pcalog)),4)
    y.hat.pcalog[which(folds == j)] <- I(as.numeric(predict(pca.logistic_fit, as.data.frame(t(jfolddata[,1])), type = "response")) > 0.5)
  }
  
  error.results[j+15,] <- c(
    mean(I(y.hat.pcalog[which(folds == j)] != t(Y)[which(folds == j)])),
    4)
  
  # elastic net (ridge regressions) ---------------------------------------------------------------------------------------------
  
  enet.fit <<- train(x = X[-which(folds == j),], y = as.factor(Y[-which(folds == j),]), 
                    method = "glmnet", trControl = cv_5, tuneGrid = ridge.grid)
  
  cat("\n", names(enet.fit))
  cat("\n", enet.fit$pred)
  y.hat.enet[which(folds == j)] <- predict(enet.fit, newdata = X[which(folds == j),])
  error.results[j+20,] <- c(
    mean(I(y.hat.enet[which(folds == j)] != t(Y+1)[which(folds == j)])),
    5)
}

y.hat.enet <- y.hat.enet - 1
#lambda.hat.AIC
#lambda.hat.BIC

#cbind(y.hat.AIC, Y[seq(1, length(LOO_set))])
y.true <- Y

finalerror <- c(
  mean(I(y.hat.AIC != t(y.true))),
  mean(I(y.hat.BIC != t(y.true))),
  mean(I(y.hat.log != t(y.true))),
  mean(I(y.hat.pcalog != t(y.true))),
  mean(I(y.hat.enet != t(y.true)))
)

finalerror <- t(as.data.frame(finalerror))

colnames(finalerror) <- c("RgPCC.AIC", "RgPCC.BIC", "logistic", "PCA.logistic", "Ridge")

g1 <- tableGrob(finalerror)

h <- grid::convertHeight(sum(g1$heights), "mm", TRUE)
w <- grid::convertHeight(sum(g1$widths), "mm", TRUE)
ggplot2::ggsave(paste(name, ".pdf"), g1, height = h, width = w + 3, units = "mm")

end_time <- Sys.time()
cat(" time :", difftime(end_time, start_time, units = time_units), time_units, "\n")



# ===================================================================================== statistical analysis

anova <- aov(cverror ~ as.factor(method), data = error.results)
tukey <- TukeyHSD(anova, which = "as.factor(method)")

# plot of errors

ggplot(error.results, aes(x=method, y=cverror)) + 
    geom_point()+
    stat_summary(aes(y = cverror,group=1), fun=mean, colour="red", geom="line",group=1)+
    stat_summary(fun.data = mean_se, geom = "errorbar", colour="red")+
    labs(title = name, caption = paste("1 = RgPCC.AIC, 2 = RgPCC.BIC, 3 = logreg, 4 = pcalogreg.90 ( #PC = ", num.pc,  "), 5 = ridge"))
ggsave(paste(name, "graph.pdf"), height = 12, width = 16, units = "cm")

# plot the AIC data

colnames(AIC.data) <- c("AIC", "lambda", "fold")
ggplot(AIC.data, aes(x = lambda, y = AIC, group = fold, color = as.factor(fold))) +
  geom_point() + 
  geom_line() +
  labs(color = "Folds")
ggsave(paste(name, "AIC.pdf"), height = 12, width = 16, units = "cm")

# plot the BIC data

colnames(BIC.data) <- c("BIC", "lambda", "fold")
ggplot(BIC.data, aes(x = lambda, y = BIC, group = fold, color = as.factor(fold))) +
  geom_point() + 
  geom_line() +
  labs(color = "Folds")
ggsave(paste(name, "BIC.pdf"), height = 12, width = 16, units = "cm")

# ======================================================================================== outputs

output <- list(finalerror, 
               lambda.hat.AIC, 
               lambda.hat.BIC, 
               paste("(", lambda.set[1], ", ", lambda.set[2], ", ..., ", lambda.set[length(lambda.set)], ")", sep = ""),
               error.results, 
               anova, 
               tukey)

names(output) <- c("finalerror", "lambda.hat.AIC.indices", "lambda.hat.BIC.indices", "lambda.set", "error.results", "anova", "pwtukey")

setwd(original.wd)

return(output)
}