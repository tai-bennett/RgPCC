# testing the new RgPCC function for simulated data

source('../functions/simulation_function.r')

# =================================================================================
# PARAMETERS
# =================================================================================

p <- 12 # number of predictor variables
my.sample_size_set <- c(100,200) # sample size of each sample

# covariance matrix generation
my.rho <- 0.8

# tolerance of Newton-Ralphson method
my.tol <- 0.1

# number of replications of the experiment
my.N <- 35 

# where output gets saved
mywd <- "results/script5"


# =================================================================================
# sparsity 1, debug
# =================================================================================

case1results <- RgPCC.lasso.simulated.exp.v4(sample_size_set = my.sample_size_set,
                          gamma_set = rbind(c(rep(2,5) 20, 10, 10, rep(2, p-8))),
						  gamma_index = "0",
                          rho = my.rho, 
                          p = p, 
                          N = my.N, 
                          lambda_set = seq(0, 60, 2), 
                          time_units = "min", 
                          tol_0 = my.tol,
                          wd = mywd,
                          save.name="script-5")


