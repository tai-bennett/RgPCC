# RgPCC_lasso on simulated data function v2
#
# can control paramters p (dimensionality) N (sample size) gamma.set
#
#
#
#
#

RgPCC.lasso.simulated.exp.v2 <- 
  function (sample_size_set, 
            gamma_set, 
            rho, 
            p,
            sample_size_factor = 5, 
            N, 
            lambda_set, 
            time_units = "sec", 
            tol_0 = 0.1,
            wd,
            save.name){
    
    setwd(wd)
    
    # start ========================================================================================= PACKAGES
    library(MASS)
    library(Rlab)
    library(matlib)
    library(ggplot2)
    library(stargazer)
    library(grid)
    library(gridExtra)
    # end =========================================================================================== PACKAGES
    
    #source("D:/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/RgPCC_lasso_experimental.R")
    #source("D:/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/RgPCC_lasso_experimental_v2.R")
    #source("D:/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/RTG_functions.R")
    
	source('RgPCC_lasso_experimental.R')
	source('RgPCC_lasso_experimental_v2.R')
    source('RTG_functions.R')
	
    
    # start ========================================================================================= INITIALIZE
    
    gamma.index.set <- seq(1, nrow(gamma_set), 1) # set of indices for gammas to be tested
    lambda.index.set <- seq(1,length(lambda_set), 1) # set of indices for the lambdas
    
    Sigma <- matrix(c(1:(p*p)), nrow = p); # Sigma covariance matrix
    for(i in 1:p) {
      for (j in 1:p) {
        Sigma[i,j] <- (rho)^abs(i - j)
      }
    }
    
    # FINAL RESULTS STORAGE
    result_size <- length(sample_size_set)*length(lambda_set)
    result_size_row <- length(sample_size_set)*7
    
    AIC.df <- data.frame(
      sample_size = rep(0, result_size),
      lambda = rep(0, result_size),
      AIC = rep(0, result_size),
      stringsAsFactors = FALSE)
    BIC.df <- data.frame(
      sample_size = rep(0, result_size),
      lambda = rep(0, result_size),
      BIC = rep(0, result_size),
      stringsAsFactors = FALSE)
    MSE.df <- data.frame(
      sample_size = rep(0, result_size),
      lambda = rep(0, result_size),
      MSE = rep(0, result_size),
      stringsAsFactors = FALSE)
    pMSE.df <- data.frame(
      sample_size = rep(0, result_size),
      lambda = rep(0, result_size),
      pMSE = rep(0, result_size),
      stringsAsFactors = FALSE)
    MSECV.df <- data.frame(
      sample_size = rep(0, result_size),
      lambda = rep(0, result_size),
      MSECV = rep(0, result_size),
      stringsAsFactors = FALSE)
    
    results.df <- data.frame(
      method = rep("NA", result_size_row),
      sample_size = rep(0, result_size_row),
      one.norm = rep(0, result_size_row),
      two.norm = rep(0, result_size_row),
      EGKL = rep(0, result_size_row),
      class.error = rep(0, result_size_row),
      gamma.size = rep(0, result_size_row),
      stringsAsFactors = FALSE)
    
    
    # end +========================================================================================== INITIALIZE
    
    
    # start ========================================================================================= THE LOOPS
    total_time_start <- Sys.time()
    for (c in gamma.index.set) { # GAMMA LOOP - c index through gamma set
      for(n in sample_size_set) { # SAMPLE SIZE LOOP - n is sample size
        
        # storage matrices
        my.metrics.logistic <- matrix(0, N, 10)
        # tuning matrices
        AIC.mat <- matrix(0, N, length(lambda_set))
        BIC.mat <- matrix(0, N, length(lambda_set))
        MSE.mat <- matrix(0, N, length(lambda_set))
        pMSE.mat <- matrix(0, N, length(lambda_set))
        MSECV.mat <- matrix(0, N, length(lambda_set))
        
        # make tuning and testing data for sample size n
        data.tune <- sim_data(n, rep(0, p), Sigma, gamma_set[c,], n + 10)
        data.test <- sim_data(sample_size_factor * n, rep(0, p), Sigma, gamma_set[c,], n + 1)
        
        all.RgPCC.results <- list()
        # all.RgPCC.results[[N]][[lambda_index]]
        for (i in c(1:N)) { # N replicates of the data
          max_iterations <- 0
          if (i == 1) {start_time <- Sys.time(); cat("|")}
          cat("=")
          if (i == N/2) {cat("|")}
          
          # make ith replicated data data
          data.train <- sim_data(n, rep(0,p), Sigma, gamma_set[c,], 4*i)
          
          # run log and pca log for data
          logistic.results <- logistic.process(data.train$X, 
                                               data.train$Y, 
                                               data.train$prob, 
                                               data.test$X, 
                                               data.test$Y, 
                                               data.test$prob)
          
          
          
          my.metrics.logistic[i,] <- c(
            logistic.results$mymetrics.test.log,
            logistic.results$class.error.test.log,
            logistic.results$gamma.size.log,
            logistic.results$mymetrics.test.logpca,
            logistic.results$class.error.test.logpca,
            logistic.results$gamma.size.pcalog
          )
          
          RgPCC.results.of.replicate <- list()
          for(lambda.index in lambda.index.set) { # LAMBDA LOOP - lambda tuning parameter
            
            # train for each lambda
            RgPCC_lasso_results <- RgPCC_lasso_experimental_v2(
              data.train$X, data.train$Y, lambda_set[lambda.index], tol_0, Sigma
            )
            RgPCC.results.of.replicate[[length(RgPCC.results.of.replicate)+1]] <- RgPCC_lasso_results
            
            # store AIC, BIC, MSE etc
            loglikelihood <- sum(data.train$Y * log(RgPCC_lasso_results$prob_hat) 
                                 + (1-data.train$Y) * log(1-RgPCC_lasso_results$prob_hat))
            
            AIC.mat[i, lambda.index] <- 2 * RgPCC_lasso_results$gamma_size - 2 * loglikelihood
            BIC.mat[i, lambda.index] <- RgPCC_lasso_results$gamma_size * log(n) - 2 * loglikelihood
            MSE.mat[i, lambda.index] <- RgPCC.predict(data.tune$X, data.tune$Y, RgPCC_lasso_results$beta_hat)$MSE
            pMSE.mat[i, lambda.index] <- get.my.metrics(RgPCC.predict(data.tune$X, data.tune$Y, RgPCC_lasso_results$beta_hat)$p.hat, data.tune$prob)[[2]]
            MSECV.mat[i, lambda.index] <- get.MSECV(5, data.train$X, data.train$Y, lambda_set[lambda.index], tol_0, Sigma)$CVerror
          }
          
          all.RgPCC.results[[length(all.RgPCC.results)+1]] <- RgPCC.results.of.replicate
          
          #print(lapply(list.of.RgPCC.results, "[[", 2))
          
          #after tuning we get metrics for this replicate
          
          if (i == N) {end_time <- Sys.time()
          cat("| time :", difftime(end_time, start_time, units = time_units), time_units, "\n")}
        } # loop of N replicates
        
        tuned.lambda <- get.my.lambda(AIC.mat, BIC.mat, MSE.mat, pMSE.mat, MSECV.mat, lambda_set)
        
        best.AIC <- lapply(all.RgPCC.results, "[[", tuned.lambda$AIC.index)
        best.BIC <- lapply(all.RgPCC.results, "[[", tuned.lambda$BIC.index)
        best.MSE <- lapply(all.RgPCC.results, "[[", tuned.lambda$MSE.index)
        best.pMSE <- lapply(all.RgPCC.results, "[[", tuned.lambda$pMSE.index)
        best.MSECV <- lapply(all.RgPCC.results, "[[", tuned.lambda$MSECV.index)
        
        
        adjust <- (match(n, sample_size_set)-1)*length(lambda_set)
        adjust2 <- (match(n, sample_size_set)-1)*7
        current.entries <- 1:length(lambda_set)+adjust
        
        
        AIC.df[current.entries, ] <- store.tuning.data(n, lambda_set, AIC.mat)
        BIC.df[current.entries, ] <- store.tuning.data(n, lambda_set, BIC.mat)
        MSE.df[current.entries, ] <- store.tuning.data(n, lambda_set, MSE.mat)
        pMSE.df[current.entries, ] <- store.tuning.data(n, lambda_set, pMSE.mat)
        MSECV.df[current.entries, ] <- store.tuning.data(n, lambda_set, MSECV.mat)
        
        results.df[1+adjust2,] <- c("log", n, round(colMeans(my.metrics.logistic)[1:5], 4))
        results.df[2+adjust2,] <- c("pcalog", n, round(colMeans(my.metrics.logistic)[6:10], 4))
        results.df[3+adjust2,] <- c("RgPCC.AIC", n, round(metrics(best.AIC, N, data.test$X, data.test$Y, data.test$prob), 4))
        results.df[4+adjust2,] <- c("RgPCC.BIC", n, round(metrics(best.BIC, N, data.test$X, data.test$Y, data.test$prob), 4))
        results.df[5+adjust2,] <- c("RgPCC.MSE", n, round(metrics(best.MSE, N, data.test$X, data.test$Y, data.test$prob), 4))
        results.df[6+adjust2,] <- c("RgPCC.pMSE", n, round(metrics(best.pMSE, N, data.test$X, data.test$Y, data.test$prob), 4))
        results.df[7+adjust2,] <- c("RgPCC.MSECV", n, round(metrics(best.MSECV, N, data.test$X, data.test$Y, data.test$prob), 4))
      }
    }
    
    # make tuning graphs
    
    tune.visuals(AIC.df, BIC.df, MSE.df, pMSE.df, MSECV.df, save.name)
    metric.visuals(results.df, save.name)
    
  }# end of function






























