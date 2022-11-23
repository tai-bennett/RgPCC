results <- readRDS(file = "../results/script6/results.rds")

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
}
# for each weighted data compute eigenbasis to compare to eigenbasis of X^T X
