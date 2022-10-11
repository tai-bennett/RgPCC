# =============================================================================
# =============================================================================
# =============================================================================
# Title : RgPCC_lasso simulated v3
# Author : Duncan Bennett
# -----------------------------------------------------------------------------
# Code Description : this code runs RgPCC_lasso_simualated_experiment_function
# v2 on a set of parameters that vary dimension, sparsity etc.
# =============================================================================
# =============================================================================
# =============================================================================

source("RTG_functions.R")

# Sigma covariance matrix
rho <- 0.8
p <- 2
Sigma <- matrix(c(1:(p*p)), nrow = p); 
for(i in 1:p) {
  for (j in 1:p) {
    Sigma[i,j] <- (rho)^abs(i - j)
  }
}

mydata <- sim_data(
  100, 
  rep(0, 2), 
  Sigma, 
  c(1,10), 
  123123123)

mydata2 <- sim_data(
  100,
  rep(0,2),
  Sigma,
  c(30,1),
  123123132
)

X <- data.frame(mydata[[1]])
Y <- data.frame(mydata[[2]])
P <- data.frame(mydata[[3]])

classData <- cbind(X,Y)
colnames(classData) <- c("X1", "X2", "Y")
regData <- cbind(X,P)
colnames(regData) <- c("X1", "X2", "P")

g <- ggplot(data = data.frame(classData), aes(x = X1, y = X2)) + geom_point() +
  geom_segment(aes(x = -3, y = -3, xend = 3, yend = 3)) + 
  geom_segment(aes(x = -2, y = 2, xend = 2, yend = -2))
print(g)

g <- ggplot(data = data.frame(classData), aes(x = X1, y = X2)) + geom_point(aes(color = as.factor(Y))) +
  geom_segment(aes(x = -3, y = -3, xend = 3, yend = 3)) + 
  geom_segment(aes(x = -2, y = 2, xend = 2, yend = -2))
print(g)

g <- ggplot(data = data.frame(regData), aes(x = X1, y = X2)) + geom_point(aes(color = P)) +
  scale_colour_gradientn(colours = terrain.colors(10)) + 
  geom_segment(aes(x = -3, y = -3, xend = 3, yend = 3)) + 
  geom_segment(aes(x = -2, y = 2, xend = 2, yend = -2))
print(g)

X <- data.frame(mydata2[[1]])
Y <- data.frame(mydata2[[2]])
P <- data.frame(mydata2[[3]])

classData <- cbind(X,Y)
colnames(classData) <- c("X1", "X2", "Y")
regData <- cbind(X,P)
colnames(regData) <- c("X1", "X2", "P")

g <- ggplot(data = data.frame(classData), aes(x = X1, y = X2)) + geom_point(aes(color = as.factor(Y))) +
  geom_segment(aes(x = -3, y = -3, xend = 3, yend = 3)) + 
  geom_segment(aes(x = -2, y = 2, xend = 2, yend = -2))
print(g)

g <- ggplot(data = data.frame(regData), aes(x = X1, y = X2)) + geom_point(aes(color = P)) +
  scale_colour_gradientn(colours = terrain.colors(10)) + 
  geom_segment(aes(x = -3, y = -3, xend = 3, yend = 3)) + 
  geom_segment(aes(x = -2, y = 2, xend = 2, yend = -2))
print(g)

