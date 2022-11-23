# testing here
source("../functions/all.r")

my_p <- 12 # number of predictor variables
my_sample_size <- 300 # sample size of each sample

# covariance matrix generation
my_rho <- 0.8

# tolerance of Newton-Ralphson method
my_tol <- 0.1

# number of replications of the experiment
my_N <- 35

my_gamma <- c(rep(2, 5), 20, 10, 10, rep(2, my_p - 8))
# my_gamma <- c(rep(0, 5), 20, 10, 10, rep(0, my_p - 8))

# where output gets saved
results <- tuner_rgpcc(
    my_sample_size,
    my_p,
    my_rho,
    my_gamma
)

results
save_name <- paste("../results/script6/results-nonsparse.rds", sep = "")
# save_name <- paste("../results/script6/resultse.rds", sep = "")
saveRDS(results, file = save_name)
