
tuner_rgpcc <- function(sample_size, p, rho, gamma, seed = 1932128, k = 4) {

    sigma <- matrix(c(1:(p*p)), nrow = p); 
    for(i in 1:p) {
      for (j in 1:p) {
        sigma[i,j] <- (rho)^abs(i - j)
      }
    }

    # make tuning and testing data for sample size n
    data_train <- sim_data(
      sample_size,
      rep(0, p),
      sigma,
      gamma,
      seed)
    data_tune <- sim_data(
      sample_size,
      rep(0, p),
      sigma,
      gamma,
      seed + 23984)
    data_test <- sim_data(
      sample_size,
      rep(0, p),
      sigma,
      gamma,
      seed + 284)

    data_dict <- list(
        "train" = data_train,
        "tune" = data_tune,
        "test" = data_test
    )

    grid_step <- 10
    grid <- list(
              "AIC" = 50 + seq(-50, 50, by = grid_step),
              "BIC" = 50 + seq(-50, 50, by = grid_step),
              "MSE" = 50 + seq(-50, 50, by = grid_step),
              "pMSE" = 50 + seq(-50, 50, by = grid_step),
              "MSECV" = 50 + seq(-50, 50, by = grid_step)
    )
    grid_size <- length(grid[["AIC"]])
    tuning_data <- list(length = 5 * k)
    l <- 1

    best_lambdas <- list()
    best_models  <- list()
    
    for (method in names(grid)) {
        print(paste("Starting method", method, sep=" "))
        grid_step <- 10
        for (i in 1:k) {
            # get metrics for each lambda in current grid
            tuning_results <- rep(0, length(grid[[method]]))
            for (j in 1:length(grid[[method]])) {
                lambda <- grid[[method]][j]
                # fit rgpcc for each lambda and store tuning parameters
                results <-
                    RgPCC_lasso_experimental_v3(
                                                data_train$X,
                                                data_train$Y,
                                                grid[[method]][j],
                                                tol = 0.1)

                # compute tuner metrics for each metric
                tuning_results[j] <- tuning_metric(data_train, results, data_tune, metric = method, lambda = lambda)
            }
            print(tuning_results)
            # store tuning results
            tuning_data[[l]] <-
                t(rbind(grid[[method]], tuning_results, rep(method, length(tuning_results))))
            l <- l + 1
            # make new grid
            best_idx <- which.min(tuning_results)
            best_lambda <- grid[[method]][best_idx]
            grid_step <- grid_step / 10
            grid[[method]] <-
                best_lambda + seq(-5 * grid_step, 5 * grid_step, by = grid_step)
            # ensure only nonnegative lambda
            grid[[method]] <- grid[[method]][!grid[[method]] < 0]

            if (i == k) {
                best_lambdas[[method]] <- best_lambda
                best_models[[method]] <-
                    RgPCC_lasso_experimental_v3(
                                                data_train$X,
                                                data_train$Y,
                                                best_lambda,
                                                tol = 0.1)
            }
        }
    }
    tuning_data <- do.call(rbind, tuning_data)
    tuning_data[,1] <- as.double(tuning_data[,1])
    tuning_data[,2] <- as.double(tuning_data[,2])
    output <- list(
        "tuning_data" = tuning_data,
        "best_lambda" = best_lambdas,
        "best_models" = best_models,
        "data_dict" = data_dict
    )
}

tuning_metric <- function(data_train, results, data_tune, metric, lambda) {
    if (metric == "AIC") {
        loglikelihood <-
            sum(data_train$Y * log(results$prob_hat)
                             + (1 - data_train$Y) * log(1 - results$prob_hat))

        output <- 2 * results$gamma_size - 2 * loglikelihood
    }
    if (metric == "BIC") {
        sample_size <- dim(data_train$X)[1]
        loglikelihood <-
            sum(data_train$Y * log(results$prob_hat)
                             + (1 - data_train$Y) * log(1 - results$prob_hat))

        output <- results$gamma_size * log(sample_size) - 2 * loglikelihood
    }
    if (metric == "MSE") {
        output <-
            RgPCC.predict(data_tune$X, data_tune$Y, results$beta_hat)$MSE
    }
    if (metric == "pMSE") {
        output <-
            get.my.metrics(
                           RgPCC.predict(
                                         data_tune$X,
                                         data_tune$Y,
                                         results$beta_hat
                           )$p.hat,
                           data_tune$prob
                           )[[2]]
    }
    if (metric == "MSECV") {
        output <-
            get.MSECValt(
                         5,
                         data_train$X,
                         data_train$Y,
                         lambda,
                         tol_0 = 0.1,
                         sigma
            )$CVerror
    }
    return(output)
}
