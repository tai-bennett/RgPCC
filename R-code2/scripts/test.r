library(ggplot2)
print(results$tuning_data)

lambda <- as.double(results$tuning_data[, 1])
score <- as.double(results$tuning_data[, 2])
method <- results$tuning_data[, 3]

mydata <- data.frame(lambda, score, method, stringsAsFactors = TRUE)

mydata_info <- mydata[mydata$method == "AIC" | mydata$method == "BIC"]

pdf("tuning_plot.pdf")
myplot <- ggplot(mydata, aes(x = lambda, y = score, color = method)) +
    geom_line() + geom_point()
print(myplot)
dev.off()
