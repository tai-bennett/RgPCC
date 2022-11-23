library(ggplot2)
nonsparse <- TRUE

if (nonsparse == TRUE) {
results <- readRDS(file = "../results/script6/results-nonsparse.rds")
} else {
results <- readRDS(file = "../results/script6/results.rds")
}

# ===========================================================================
# finding angles between spaces
# ===========================================================================
# extract things from the results object
method_set <- c("AIC", "BIC", "MSE", "pMSE", "MSECV")
# for each tuning method, get weights and computed weighted data
for (method in method_set) {
    data <- results$data_dict$train
    original_svd <- svd(data$X)
    weighted_svd <- svd(results$best_models[[method]]$W_final %*% data$X)
    # get 90% of the variance for each subspace
    Do <- original_svd$d
    So <- sum(Do)
    Dw <- weighted_svd$d
    Sw <- sum(Dw)
    partial_So <- 0
    partial_Sw <- 0
    i_So <- 0
    i_Sw <- 0
    while (partial_So < 0.9) {
        i_So <- i_So + 1
        partial_So <- partial_So + Do[i_So]/So
    }
    while (partial_Sw < 0.9) {
        i_Sw <- i_Sw + 1
    partial_Sw <- partial_Sw + Dw[i_Sw]/Sw
    }

    Vo <- original_svd$v[,1:i_So]
    Vw <- weighted_svd$v[,1:i_Sw]

    S <- svd(t(Vo) %*% Vw)$d
    print(paste(
                "method = ", method,
                " with ",
                " lambda = ", results$best_lambda[[method]],
                " and ",
                "nonsparse = ", nonsparse,
                sep = ""))
    print(S)
}
# for each weighted data compute eigenbasis to compare to eigenbasis of X^T X

# ===========================================================================
# ploting tuning data
# ===========================================================================
df <- data.frame(results$tuning_data)
df[, 1] <- as.double(df[, 1])
df[, 2] <- as.double(df[, 2])
colnames(df) <- c("lambda", "metric", "method")
for (i in 1:2) {
    if (i == 1) {
        df2 <- df[df$method == "AIC" | df$method == "BIC", ]
    }
    if (i == 2) {
        df2 <- df[df$method == "MSE" | df$method == "pMSE" | df$method == "MSECV", ]
    }
    plt <- ggplot(data = df2, aes(x = lambda, y = metric, color = method)) +
        geom_point() + geom_line()
    if (nonsparse == TRUE) {
        pdf(paste("../results/script7/tuning_results-nonsparse", i, ".pdf", sep = ""))
    } else {
        pdf(paste("../results/script7/tuning_results", i, ".pdf", sep = ""))
    }
    print(plt)
    dev.off()
}
