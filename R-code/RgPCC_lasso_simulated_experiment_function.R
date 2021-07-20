# RgPCC_lasso on simulated data function
#
# can control paramters p (dimensionality) N (sample size) gamma.set
#
#
#
#
#

RgPCC.lasso.simulated.exp <- 
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
    
    source('RgPCC_lasso_experimental.R')
	source('RgPCC_lasso_experimental_v2.R')
    source('RTG_functions.R')
    
    
    # start ========================================================================================= F-RESULTS DATAFRAME
    result_size <- 1000
    results <- data.frame(
      gamma_true = rep(0, result_size),
      sample_size = rep(0, result_size),
      lambda = rep(0, result_size),
      MSE_prob = rep(0, result_size),
      MSE_train = rep(0, result_size),
      MSE_gamma_hat = rep(0, result_size),
      MSE_beta_hat = rep(1, result_size),
      MSE_prob_tune = rep(0, result_size),
      MSE_tune = rep(0, result_size),
      AIC = rep(0,result_size),
      BIC = rep(0,result_size),
      loglike = rep(0,result_size),
      MSECV = rep(0,result_size),
      MSE_prob_test = rep(0, result_size),
      MSE_test = rep(0, result_size),
      mean_gamma_size = rep(0, result_size),
      max_iter = rep(0,result_size),
      mean_iter = rep(0,result_size),
      stringsAsFactors = FALSE)
    
    
    results_misc <- data.frame(
      gamma_true = rep(0, result_size),
      sample_size = rep(0, result_size),
      Bayes_test_error = rep(0, result_size),
      MSE_train_log = rep(0, result_size),
      MSE_test_log = rep(0, result_size),
      MSE_prob_train_log = rep(0, result_size),
      MSE_prob_test_log = rep(0, result_size),
      MSE_train_pca_log = rep(0, result_size),
      MSE_test_pca_log = rep(0, result_size),
      MSE_prob_train_pca_log = rep(0, result_size),
      MSE_prob_test_pca_log = rep(0, result_size)
    )
    
    
    
    
    # end =========================================================================================== F-RESULTS DATAFRAME
    
    
    
    
    # start ========================================================================================= DEFINE PARAMETER
    # ------------------------------------- misc parameters
    # ------------------------------------- data generation parameters

    # set of true gammas
    #gamma_set <- rbind(
    #  c(25, rep(0, p-1)),
    #  c(20, 10, 10, rep(0, p-3)),
    #  c(15, 10, 5, 5, 3, rep(0, p-5)),
    #  c(20, 0, 0.4, 0.8, 0.8, 1.2, 0.4, rep(0, p-7)),
    #  c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, rep(0, p-12)),
    #  c(0, 0, 0.4, 0.8, 0.8, 1.2, 0.4, 0, 0, 0, 0, 20, rep(0, p-12))
    #)
    
    # set of indices for gammas to be tested
    gamma_index_set <- seq(1, nrow(gamma_set), 1)
    
    # covariance matrix generation
    #rho <- 0.8
    
    # Sigma covariance matrix
    Sigma <- matrix(c(1:(p*p)), nrow = p);
    for(i in 1:p) {
      for (j in 1:p) {
        Sigma[i,j] <- (rho)^abs(i - j)
      }
    }
    
    # create set of true betas
    v_true <- eigen(Sigma)$vector
    d_true <- diag(eigen(Sigma)$values)
    # true beta from simulated data beta = U%*%gamma (dependent on each simulation)
    beta_set <- matrix(0, nrow(gamma_set), p)
    for (r in seq(nrow(gamma_set))) {
      beta_set[r, 1:p] <- t(v_true) %*% inv(sqrt(d_true)) %*% gamma_set[r, 1:p]
    }
    
    # ------------------------------------- parameters to vary over each test
    # tolerance of Newton-Ralphson method
    #tol_0 <- 0.1
    
    # tuning parameters and sample sizes
    #lambda_set <- c(50)
    
    #sample_size_factor <- 5 # how much larger the testing set is than the training set
    
    #N <- 100 # N is number of samples
    
    count <- 0 # counter for the progress of total iterations
    
    # end =========================================================================================== DEFINE PARAMETERS
    
    
    
    
    # start ========================================================================================= THE LOOPS
    total_time_start <- Sys.time()
    for (c in gamma_index_set) { # GAMMA LOOP - c index through gamma set
      for(n in sample_size_set) { # SAMPLE SIZE LOOP - n is sample size
        for(lambda in lambda_set) { # LAMBDA LOOP - lambda tuning parameter
          ##------------------------------------------------------------------------- Progress Update
          count <- count + 1
          cat("\n", "(gamma, size, lambda) = ", "(", c,",", n,",", lambda, ")")
          cat("\n", count) 
          cat(" of", length(gamma_index_set) * length(sample_size_set) * length(lambda_set), "\n")
          
          ##------------------------------------------------------------------------- Empty Matrices for Gamma, Beta, Prob
          
          ###.......................................... predicted
          gamma_hat <- matrix(0, nrow = N, ncol = p) # each row is a predicted gamma
          beta_hat <- matrix(0, N, p) # each row is a predicted beta
          prob_hat <- matrix(0, N, n) # each row is a predicted probabilities of class 1 (training data)
          prob_hat_test <- matrix(0, 1, sample_size_factor * n) # predicted probabilities of class 1 (testing data)
          prob_hat_tune <- matrix(0, 1, n) 
          
          ###.......................................... error
          tol_data <- matrix(0, N, 1) # tolerance (percent change of betas) when RgPCC_lasso was aborted
          MSE_prob <- matrix(0, N, 1) # MSE error of predicted probabilities (training data)
          MSE_prob_test <- matrix(0, N, 1) # MSE error of predicted probabilities (testing data)
          MSE_prob_tune <- matrix(0, N, 1)
          MSE_gamma_hat <- matrix(0, N, 1) # MSE of the gamma hats
          MSE_beta_hat <- matrix(0, N, 1) # MSE of the beta hats
          MSECV <- matrix(0, N, 1)
          MSECV_se <- matrix(0, N, 1)
          error_test <- matrix(0, N, 1)
          error_tune <- matrix(0, N, 1)
          error_train <- matrix(0, N, 1)
          error_prob_test_log <- matrix(0, N, 1)
          error_prob_train_log <- matrix(0, N, 1)
          error_test_log <- matrix(0, N, 1)
          error_train_log <- matrix(0, N, 1)
          error_prob_test_pca_log <- matrix(0, N, 1)
          error_prob_train_pca_log <- matrix(0, N, 1)
          error_test_pca_log <- matrix(0, N, 1)
          error_train_pca_log <- matrix(0, N, 1)
          bayes_test_error_set <- matrix(0, N, 1)
          
          my_metrics_log <- matrix(0, N, 3)
          
          my_metrics_pca_log <- matrix(0, N, 3)
          
          
          ###.......................................... tuning
          AIC <- matrix(0, N, 1)
          BIC <- matrix(0, N, 1)
          loglike <- matrix(0, N, 1)
          
          
          ###.......................................... truth
          gamma_true <- matrix(0, N, p) # true gamma (based on simulated data)
          prob_true <- matrix(0, N, n)  # true prob based on each simulation (logistic(X%*%beta_true)
          
          
          ###.......................................... misc
          
          gamma_sizes <- matrix(0, N, 1) # store sizes of gamma_hat
          iteration_set <- matrix(0, N, 1) # iteration set
          
          ##------------------------------------------------------------------------- Tuning Data
          data_tune <- sim_data(
            n, rep(0, p), Sigma, gamma_set[c,], n + 10
          )
          prob_true_tune <- data_tune[[3]]
          
          
          ##------------------------------------------------------------------------- Testing Data
          data_test <- sim_data(
            sample_size_factor * n, rep(0, p), Sigma, gamma_set[c,], n + 1
          )
          
          ###.......................................... testing data info
          Y_hat_bayes <- I(data_test[[1]] %*% beta_set[c, ] > 0.5)
          prob_true_test <- data_test[[3]]
          
          
          ## ------------------------------------------------------------------------- *Samples Loop*
          
          for (i in c(1:N)) { # SAMPLES LOOP - i indexes through the N replicates
            max_iterations <- 0
            ###.......................................... progress of samples
            if (i == 1) {start_time <- Sys.time(); cat("|")}
            if (0 == i%%2) {cat("=")}
            if (0 == i%%50) {cat("|")}
            
            ###.......................................... create training data
            data_train <- sim_data(n, rep(0,p), Sigma, gamma_set[c,], 4*i)
            
            ###.......................................... partition data for cross validation
            
            MSE_parts <- c(0,0,0,0,0)
            prob_hat_cv <- matrix(0, 1, n)
            for (k in c(1:5)) {
              cv <- c((n/5 * (k - 1) + 1):(n/5 * k))
              X_cv <- data_train[[1]][-cv, ]
              Y_cv <- data_train[[2]][-cv, ]
              RgPCC_lasso_results_cv <- RgPCC_lasso_experimental(
                X_cv, Y_cv, lambda, tol_0, Sigma
              )
              
              # estimate prob for test set
              for (j in cv) {
                prob_hat_cv[, j] <- mylogistic(data_train$X[j, ],
                                               as.matrix(RgPCC_lasso_results_cv[[2]]))
              }
              
            }
            
            for (k in c(1:5)){
              cv <- c((n/5 * (k - 1) + 1):(n/5 * k))
              MSE_parts[k] <- (5/n)*sum(I(prob_hat_cv[, cv] > 0.5) != t(data_train[[2]][cv,]))
            }
            MSECV_se <- sd(MSE_parts)/sqrt(5)
            #SSE[k] <- sum(I(prob_hat_cv[ ,cv] > 0.5) != t(data_train[[2]])[, cv])
            MSECV[i,1] <- (1/n)*sum(I(prob_hat_cv > 0.5) != t(data_train[[2]]))
            #MSECV[i,1] <- (1/n)*sum(SSE)
            
            
            ###.......................................... training data info
            Y_best <- I(data_train[[1]] %*% beta_set[c, ] > 0.5)
            prob_true <- data_train[[3]]
            
            ###.......................................... training RgPCC_lasso and results
            RgPCCtime.start <- Sys.time()
            RgPCC_lasso_results <- RgPCC_lasso_experimental(
              data_train[[1]], data_train[[2]], lambda, tol_0, Sigma
            )
            RgPCCtime.end <- Sys.time()
            #cat("RgPCC time :", difftime(RgPCCtime.end, RgPCCtime.start, units = time_units), time_units, "\n")
            
            #if (RgPCC_lasso_results[[7]]) {next}
            
            #### ____________________ storing RgPCC training prediction
            gamma_hat[i, 1:p] <- as.matrix(RgPCC_lasso_results[[1]])
            beta_hat[i, 1:p] <- as.matrix(RgPCC_lasso_results[[2]])
            prob_hat[i, 1:n] <- as.matrix(RgPCC_lasso_results[[3]])
            
            #### ____________________ storing RgPCC training errors
            MSE_prob[i, 1] <- mean((prob_true - prob_hat[i, 1:n])^2)
            MSE_gamma_hat[i, 1] <- sum((gamma_set[c, ] - gamma_hat[i,1:p])^2)
            MSE_beta_hat[i, 1] <- sum(((beta_set[c, ] - beta_hat[i,1:p])^2))
            y_hat_train <- I(prob_hat > 0.5)
            error_train[i, 1] <- mean(I(y_hat_train[i, ] != t(data_train[[2]])))
            
            #### ____________________ storing RgPCC training misc
            gamma_sizes[i, 1] <- colSums(as.matrix(gamma_hat[i, 1:p]) != 0)
            tol_data[i, 1] <- as.matrix(RgPCC_lasso_results[[4]])
            
            #### ____________________ storing RgPCC training tuning
            
            loglikelihood <- sum(data_train[[2]] * log(prob_hat[i, 1:n]) + (1-data_train[[2]]) * log(1-prob_hat[i, 1:n]))
            #        sum yi * log(Xi beta) + (1 - yi)log(Xi beta)
            AIC[i, 1] <- 2 * gamma_sizes[i, 1] - 2 * loglikelihood
            BIC[i, 1] <- gamma_sizes[i, 1] * log(n) - 2 * loglikelihood
            loglike[i,1] <- loglikelihood
            
            ###.......................................... tuning RgPCC_lasso and results
            for (j in seq_len(nrow(data_tune$X))) {
              prob_hat_tune[, j] <- mylogistic(data_tune$X[j, ],
                                               as.matrix(RgPCC_lasso_results[[2]]))
            }
            y_hat_tune <- I(prob_hat_tune > 0.5)
            error_tune[i, 1] <- mean(I(y_hat_tune != t(data_tune[[2]])))
            MSE_prob_tune[i, 1] <- mean((prob_true_tune - t(prob_hat_tune))^2)
            
            
            
            ###.......................................... testing RgPCC_lasso and results
            for (j in seq_len(nrow(data_test$X))) {
              prob_hat_test[, j] <- mylogistic(data_test$X[j, ],
                                               as.matrix(RgPCC_lasso_results[[2]]))
            }
            y_hat_test <- I(prob_hat_test > 0.5)
            error_test[i, 1] <- mean(I(y_hat_test != t(data_test[[2]])))
            MSE_prob_test[i, 1] <- mean((prob_true_test - t(prob_hat_test))^2)
            
            # save bayes error
            Bayes_error <- mean(I(prob_true_test > 0.5) != data_test[[2]])
            bayes_test_error_set[i, 1] <- Bayes_error
            
            if (RgPCC_lasso_results[[5]] == TRUE) {max_iterations <- max_iterations + 1}
            #if (count == TRUE | RgPCC_lasso_results[[5]] == TRUE) {
            #  cat("\n test", i, "sample size", n, "lambda", lambda)
            #}
            iteration_set[i, 1] <- RgPCC_lasso_results[[6]]
            
            ###.......................................... training and testing logistic regression
            log_fit_data.frame <- as.data.frame(cbind(data_train[[2]], data_train[[1]]))
            log_fit_test_data.frame <- as.data.frame(cbind(data_test[[2]], data_test[[1]]))
            
            logistictime.start <- Sys.time()
            logistic_fit <- glm(V1 ~ . -1 , data=log_fit_data.frame, family=binomial(link="logit"))
            logistictime.end <- Sys.time()
            #cat("Logistic time :", difftime(logistictime.end, logistictime.start, units = time_units), time_units, "\n")
            
            p_train <- predict(logistic_fit, log_fit_data.frame[2:(p+1)], type = "response")
            p_test <- predict(logistic_fit, log_fit_test_data.frame[2:(p+1)], type = "response")
            
            loglikelihood <- sum(data_train[[2]] * log(p_train) + (1-data_train[[2]]) * log(1-p_train))
            
            y_train <- I(p_train > 0.5)
            y_test <- I(p_test > 0.5)
            
            error_prob_train_log[i, 1] <- mean(I(y_train != t(data_train[[2]])))
            error_prob_test_log[i, 1] <- mean(I(y_test != t(data_test[[2]])))
            
            error_train_log[i, 1] <- mean((prob_true - p_train)^2)
            error_test_log[i, 1] <- mean((prob_true_test - p_test)^2)
            
            
            # my metrics for logistic reg
            
            my_metrics_log[i,] <- get.my.metrics(p_test, prob_true_test)
            
            
            ###.......................................... training and testing PCA logistic regression
            
            pca <- princomp(data_train[[1]], cor=F)
            
            cumsum.var <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
            num.pc <- min(which(cumsum.var > 0.90))
            #cat("\n", "num.pc = ", num.pc)
            pca.X <- pca$scores[,c(1:num.pc)]
            pca.YX <- as.data.frame(cbind(data_train[[2]], pca.X))
            
            pca.logistic_fit <- glm(V1 ~ . -1 , data=pca.YX, family=binomial(link="logit"))
            
            pca.X.test <- as.data.frame(data_test[[1]] %*% pca$loadings)
            pca.X.test <- pca.X.test[,c(1:num.pc)]
            
            p_train <- predict(pca.logistic_fit, as.data.frame(pca.X), type = "response")
            p_test <- predict(pca.logistic_fit, as.data.frame(pca.X.test), type = "response")
            
            y_train <- I(p_train > 0.5)
            y_test <- I(p_test > 0.5)
            
            error_prob_train_pca_log[i, 1] <- mean(I(y_train != t(data_train[[2]])))
            error_prob_test_pca_log[i, 1] <- mean(I(y_test != t(data_test[[2]])))
            
            error_train_pca_log[i, 1] <- mean((prob_true - p_train)^2)
            error_test_pca_log[i, 1] <- mean((prob_true_test - p_test)^2)
            
            my_metrics_pca_log[i,] <- get.my.metrics(p_test, prob_true_test)
            

            # INFO TO PRINT DURING INDIVIDUAL CASES
            #cat("\n", "beta hat is ", beta_hat[i, 1:p])
            #cat("\n", "gamma hat is ", gamma_hat[i, 1:p])
            #cat("\n", "MSE of gamma hat is ", MSE_gamma_hat[i, 1])
            #cat("\n", "log_beta is ", logistic_fit$coefficients)
            if (i == N) {end_time <- Sys.time()
            cat(" time :", difftime(end_time, start_time, units = time_units), time_units, "\n")
            }
          } # end SAMPLES LOOP - i indexes through the N samples
          
          ##------------------------------------------------------------------------- Storing Final Results
          
          # find the smallest unused row in the results dataframe to put new data in
          
          
          # ********** print for now
          
          cat(mean(MSECV_se))
          
          # ***********
          iteration_set[iteration_set == 50] <- 0
          k <- min(which(results$gamma_true == 0))
          results[k, ] <-
            c(c, n, lambda, 
              mean(MSE_prob), 
              mean(error_train),
              mean(MSE_gamma_hat), 
              mean(MSE_beta_hat),
              mean(MSE_prob_tune),
              mean(error_tune),
              mean(AIC),
              mean(BIC), 
              mean(loglike),
              mean(MSECV),
              mean(MSE_prob_test),
              mean(error_test), 
              mean(gamma_sizes),
              max_iterations,
              mean(iteration_set))
          
        } # end LAMBDA LOOP - lambda tuning parameter
        K <- min(which(results_misc$MSE_prob_test_log == 0))
        results_misc[K, ] <-
          c(c, n,
            mean(Bayes_error), 
            mean(error_prob_train_log), 
            mean(error_prob_test_log),
            mean(error_train_log),
            mean(error_test_log),
            mean(error_prob_train_pca_log), 
            mean(error_prob_test_pca_log),
            mean(error_train_pca_log),
            mean(error_test_pca_log)
          )
        
      } # end SAMPLE SIZE LOOP - n is sample size
    
    
    # end =========================================================================================== THE LOOPS
    
    
    
    
    # start ========================================================================================= RESULTS PRESENTATION
    
    ###################
    # overall results #
    ###################
    
    final_results <- round(results[1:k,],5); final_results
    
    g0 <- tableGrob(final_results)
    
    h <- grid::convertHeight(sum(g0$heights), "mm", TRUE)
    w <- grid::convertHeight(sum(g0$widths), "mm", TRUE)
    
    width_factor <- 1.7
    height_factor <- 0.5
    
    ggplot2::ggsave(paste("(", save.name,",", p, ")_final_results.pdf", sep = ""), g0,  height = h, width = w+3, units = "mm")
    
    final_results_testing <- results[1:k,]; final_results_testing
    
    #colMeans(my.metrics.log)
    
    
    # ------------------------ finding indices for min tuning metrics
    
    ###################
    # tuning  results #
    ###################
    
    final_results_tuning <- round(results[1:k,c(1,2,3, 8:15)],5); final_results_tuning
    
    # find minimized tuning
    
    minAICindex <- rep(0, length(sample_size_set))
    minBICindex <- rep(0, length(sample_size_set))
    minMSEindex <- rep(0, length(sample_size_set))
    minMSEprobindex <- rep(0, length(sample_size_set))
    minMSECVindex <- rep(0, length(sample_size_set))
    
    # find cell indices for table
    
    for(i in 1:length(sample_size_set)) {
      minAICindex[i] <- which.min(final_results_tuning[final_results_tuning$sample_size == sample_size_set[i],]$AIC) + (i - 1)*length(lambda_set)
      minBICindex[i] <- which.min(final_results_tuning[final_results_tuning$sample_size == sample_size_set[i],]$BIC) + (i - 1)*length(lambda_set)
      minMSEindex[i] <- which.min(final_results_tuning[final_results_tuning$sample_size == sample_size_set[i],]$MSE_tune) + (i - 1)*length(lambda_set)
      minMSEprobindex[i] <- which.min(final_results_tuning[final_results_tuning$sample_size == sample_size_set[i],]$MSE_prob_tune) + (i - 1)*length(lambda_set)
      minMSECVindex[i] <- which.min(final_results_tuning[final_results_tuning$sample_size == sample_size_set[i],]$MSECV) + (i - 1)*length(lambda_set)
    }
    
    mincolor <- c("pink1", "darkseagreen1", "darkgoldenrod1", "tomato1")
    
    g <- tableGrob(final_results_tuning)
    find_cell <- function(table, row, col, name="core-fg"){
      l <- table$layout
      which(l$t==row & l$l==col & l$name==name)
    }
    
    # ------------------------- color the min cells for each metric
    for (i in 1:length(sample_size_set)){
      indprob <- find_cell(g, minMSEprobindex[i] + 1, 5, "core-bg")
      indproblambda <- find_cell(g,  minMSEprobindex[i] + 1, 4, "core-bg")
      g$grobs[indprob][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      g$grobs[indproblambda][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      
      indMSE <- find_cell(g, minMSEindex[i] + 1, 6, "core-bg")
      indMSElambda <- find_cell(g, minMSEindex[i] + 1, 4, "core-bg")
      g$grobs[indMSE][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      g$grobs[indMSElambda][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      
      indAIC <- find_cell(g, minAICindex[i] + 1, 7, "core-bg")
      indAIClambda <- find_cell(g, minAICindex[i] + 1, 4, "core-bg")
      g$grobs[indAIC][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      g$grobs[indAIClambda][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      
      indBIC <- find_cell(g, minBICindex[i] + 1, 8, "core-bg")
      indBIClambda <- find_cell(g, minBICindex[i] + 1, 4, "core-bg")
      g$grobs[indBIC][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      g$grobs[indBIClambda][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      
      indMSECV <- find_cell(g, minMSECVindex[i] + 1, 10, "core-bg")
      indMSECVlambda <- find_cell(g, minMSECVindex[i] + 1, 4, "core-bg")
      g$grobs[indMSECV][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      g$grobs[indMSECVlambda][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
    }
    # ------------------------- #
    
    # draw the colored tabled
    h <- grid::convertHeight(sum(g$heights), "mm", TRUE)
    w <- grid::convertHeight(sum(g$widths), "mm", TRUE)
    ggplot2::ggsave(paste("(", save.name,",", p, ")_final_results_tuning.pdf", sep = ""), g, height = h, width = w + 3, units = "mm")
    
    g1 <- tableGrob(final_results[c(minMSEprobindex, minMSEindex, minAICindex, minBICindex, minMSECVindex), ])
    
    for (i in 1:length(sample_size_set)){
      g1$grobs[find_cell(g1, i+1, 9, "core-bg")][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      g1$grobs[find_cell(g1, i+1+1*length(sample_size_set), 10, "core-bg")][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      g1$grobs[find_cell(g1, i+1+2*length(sample_size_set), 11, "core-bg")][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      g1$grobs[find_cell(g1, i+1+3*length(sample_size_set), 12, "core-bg")][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      g1$grobs[find_cell(g1, i+1+4*length(sample_size_set), 14, "core-bg")][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
    }
    
    h <- grid::convertHeight(sum(g1$heights), "mm", TRUE)
    w <- grid::convertHeight(sum(g1$widths), "mm", TRUE)
    ggplot2::ggsave(paste("(", save.name,",", p, ")_final_results_tuning_optimal.pdf", sep = ""), g1, height = h, width = w + 3, units = "mm")
    
    
    ###################
    # testing results #
    ###################
    
    
    g2 <- tableGrob(final_results[c(minMSEprobindex, minMSEindex, minAICindex, minBICindex, minMSECVindex), c(1,2,3, 6:16)])
    
    for (i in 1:length(sample_size_set)){
      g2$grobs[find_cell(g2, i+1, 7, "core-bg")][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      g2$grobs[find_cell(g2, i+1+1*length(sample_size_set), 8, "core-bg")][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      g2$grobs[find_cell(g2, i+1+2*length(sample_size_set), 9, "core-bg")][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      g2$grobs[find_cell(g2, i+1+3*length(sample_size_set), 10, "core-bg")][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
      g2$grobs[find_cell(g2, i+1+4*length(sample_size_set), 12, "core-bg")][[1]][["gp"]] <- gpar(fill=mincolor[i], col = mincolor[i], lwd=5)
    }
    
    h <- grid::convertHeight(sum(g2$heights), "mm", TRUE)
    w <- grid::convertHeight(sum(g2$widths), "mm", TRUE)
    ggplot2::ggsave(paste("(", save.name,",", p, ")_final_results_tuning_optimaltrim.pdf", sep = ""), g2, height = h, width = w+3, units = "mm")
    
    
    # ------------------------ #
    
    ################
    # misc results #
    ################
    
    final_results_misc <- round(results_misc[1:K,], 5); final_results_misc
    
    g3 <- tableGrob(final_results_misc)
    h <- grid::convertHeight(sum(g3$heights), "mm", TRUE)
    w <- grid::convertHeight(sum(g3$widths), "mm", TRUE)
    ggplot2::ggsave(paste("(", save.name,",", p, ")_final_results_misc.pdf", sep = ""), g3, height = h, width = w+3, units = "mm")
    
    ##------------------------------------------------------------------------- Print datasets
    
    ##------------------------------------------------------------------------- Graphical Info
    
    pdf(paste("(", save.name,",", p, ")_graphtuningMSEprob.pdf", sep = ""))
    mycaption <- "Dashed (dotted) line shows testing error of probability vector for (pca) logistic regression. Solid line for RgPCC."
    myplot <- ggplot(data = final_results, aes(x = lambda, y = MSE_prob_test, group = interaction(sample_size, gamma_true))) +
      geom_line(aes(color = factor(sample_size))) +
      geom_point(aes(color = factor(sample_size))) +
      geom_hline(data = final_results_misc, linetype = "longdash", 
                 aes(yintercept = MSE_prob_test_log, color = factor(sample_size)))+
      geom_hline(data = final_results_misc, linetype = "dotted", 
                 aes(yintercept = MSE_prob_test_pca_log, color = factor(sample_size)))+
      labs(title = "Parameter Tuning MSE", subtitle = paste("(", save.name,",", p, ")"), 
           caption = mycaption, color = "Sample Size", y = "prob test error")
    print(myplot)
    dev.off()
    
    pdf(paste("(", save.name,",", p, ")_graphtuningAIC.pdf", sep = ""))
    myplotAIC <- ggplot(data = final_results, aes(x = lambda, y = AIC, group = sample_size)) +
      geom_line(aes(color = factor(sample_size))) +
      geom_point(aes(color = factor(sample_size))) +
      labs(title = "Parameter Tuning AIC", subtitle = paste("(", save.name,",", p, ")"), color = "Sample Size", y = "AIC")
    
    print(myplotAIC)
    dev.off()
    
    pdf(paste("(", save.name,",", p, ")_graphtuningBIC.pdf", sep = ""))
    myplotBIC <- ggplot(data = final_results, aes(x = lambda, y = BIC, group = sample_size)) +
      geom_line(aes(color = factor(sample_size))) +
      geom_point(aes(color = factor(sample_size))) +
      labs(title = "Parameter Tuning BIC", subtitle = paste("(", save.name,",", p, ")"), color = "Sample Size", y = "BIC")
    
    print(myplotBIC)
    dev.off()
    
    pdf(paste("(", save.name,",", p, ")_graphtuningloglike.pdf", sep = ""))
    myplotloglike <- ggplot(data = final_results, aes(x = lambda, y = loglike, group = sample_size)) +
      geom_line(aes(color = factor(sample_size))) +
      geom_point(aes(color = factor(sample_size))) +
      labs(title = "Parameter Tuning loglike", subtitle = paste("(", save.name,",", p, ")"), color = "Sample Size", y = "loglike")
    
    print(myplotloglike)
    dev.off()
    
    pdf(paste("(", save.name,",", p, ")_graphtuningMSECV.pdf", sep = ""))
    myplotAIC <- ggplot(data = final_results, aes(x = lambda, y = MSECV, group = sample_size)) +
      geom_line(aes(color = factor(sample_size))) +
      geom_point(aes(color = factor(sample_size))) +
      labs(title = "Parameter Tuning MSECV", subtitle = paste("(", save.name,",", p, ")"), color = "Sample Size", y = "MSECV")
    
    print(myplotAIC)
    dev.off()
    } # end GAMMA LOOP - c index through gamma set
    total_time_end <- Sys.time()
    cat("total time :", difftime(total_time_end, total_time_start, units = "min"), "min", "\n")
    
    # end =========================================================================================== RESULTS PRESENTATION
}



