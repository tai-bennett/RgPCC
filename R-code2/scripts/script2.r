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









# testing the new RgPCC function for simulated data

source('../functions/simulation_function.r')

# where output gets saved
mywd <- "results/script2"










casedebug = F 

#lead cases
case1 = F # sparsity 1 x
case2 = F # sparsity 1 x
case3 = F # sparsity 3 x
case4 = F # sparsity 5 x

#nonlead cases
case5 = F # sparsity 1 x
case6 = F # sparsity 1 x
case7 = F # sparsity 3 x
case8 = T # sparsity 5 x



# =================================================================================
# PARAMETERS
# =================================================================================

p <- 12 # number of predictor variables
my.sample_size_set <- c(50,100) # sample size of each sample

# covariance matrix generation
my.rho <- 0.8

# tolerance of Newton-Ralphson method
my.tol <- 0.1

# number of replications of the experiment
my.N <- 5



# =================================================================================
# =================================================================================
# SPARSITY IN LEADING PRINCIPAL COMPONENTS
# =================================================================================
# =================================================================================

# =================================================================================
# sparsity 1, N = 100
# =================================================================================

if (case1 == TRUE){
case1results <- RgPCC.lasso.simulated.exp.v4(sample_size_set = my.sample_size_set[1],
                          gamma_set = rbind(c(25, rep(0, p-1))),
						  gamma_index = "1",
                          rho = my.rho, 
                          p = p, 
                          N = my.N, 
                          lambda_set = seq(0, 60, 2), 
                          time_units = "sec", 
                          tol_0 = my.tol,
                          wd = mywd,
                          save.name="100-lead-new-algo")						  
}

# =================================================================================
# sparsity 1, N = 200
# =================================================================================

if (case2 == TRUE){
case2results <- RgPCC.lasso.simulated.exp.v4(sample_size_set = my.sample_size_set[2],
                          gamma_set = rbind(c(25, rep(0, p-1))),
						  gamma_index = "1",
                          rho = my.rho, 
                          p = p, 
                          N = my.N, 
                          lambda_set = seq(0, 150, 5), 
                          time_units = "sec", 
                          tol_0 = my.tol,
                          wd = mywd,
                          save.name="200-lead-new-algo")
}

# =================================================================================
# sparsity 3, N = {100, 200}
# =================================================================================

if (case3 == TRUE){
RgPCC.lasso.simulated.exp.v4(sample_size_set = my.sample_size_set,
                          gamma_set = rbind(c(25, rep(0, p-1))),
						  gamma_index = "2",
                          rho = my.rho, 
                          p = p, 
                          N = my.N, 
                          lambda_set = seq(0, 40, 2), 
                          time_units = "sec", 
                          tol_0 = my.tol,
                          wd = mywd,
                          save.name="lead-new-algo")
}

# =================================================================================
# sparsity 5, N = {100,200}
# =================================================================================

if (case4 == TRUE){
RgPCC.lasso.simulated.exp.v4(sample_size_set = my.sample_size_set,
                          gamma_set = rbind(c(15, 10, 5, 5, 3, rep(0, p-5))),
						  gamma_index = "3",
                          rho = my.rho, 
                          p = p, 
                          N = my.N, 
                          lambda_set = seq(0, 60, 2), 
                          time_units = "sec", 
                          tol_0 = my.tol,
                          wd = mywd,
                          save.name="leadnew-new-algo")
}



































# =================================================================================
# =================================================================================
# SPARSITY IN NONLEADING PRINCIPAL COMPONENTS
# =================================================================================
# =================================================================================

# =================================================================================
# sparsity 1, N = 100
# =================================================================================

if (case5 == TRUE){
  RgPCC.lasso.simulated.exp.v4(sample_size_set = my.sample_size_set[1],
                            gamma_set = rbind(c(0, 25, rep(0, p-2))),
							gamma_index = "1'",
                            rho = my.rho, 
                            p = p, 
                            N = my.N, 
                            lambda_set = seq(0, 60, 2), 
                            time_units = "sec", 
                            tol_0 = my.tol,
                            wd = mywd,
                            save.name="nonlead-alt-new-algo")
}


# =================================================================================
# sparsity 1, N = 100, p = 12
# =================================================================================

if (case6 == TRUE){
  RgPCC.lasso.simulated.exp.v4(sample_size_set = my.sample_size_set[2],
                            gamma_set = rbind(c(0, 0, 0, 0, 0, 25, rep(0, p-6))),
							gamma_index = "1'",
                            rho = my.rho, 
                            p = p, 
                            N = my.N, 
                            lambda_set = seq(0, 100, 10), 
                            time_units = "sec", 
                            tol_0 = my.tol,
                            wd = mywd,
                            save.name="nonlead-new-algo")
}


# =================================================================================
# sparsity 3, N = 100, p = 12 (1)
# c(20, 10, 10, rep(0, p-3)),
# c(15, 10, 5, 5, 3, rep(0, p-5))
# =================================================================================

if (case7 == TRUE){
  RgPCC.lasso.simulated.exp.v4(sample_size_set = my.sample_size_set,
                            gamma_set = rbind(c(0, 0, 0, 0, 0, 20, 10, 10, rep(0, p-8))),
							gamma_index = "2'",
                            rho = my.rho, 
                            p = p, 
                            N = my.N, 
                            lambda_set = seq(0, 20, 1), 
                            time_units = "sec", 
                            tol_0 = my.tol,
                            wd = mywd,
                            save.name="nonlead-new-algo")
}


# =================================================================================
# sparsity 5, N = 100, p = 12 (1)
# =================================================================================

if (case8 == TRUE){
  RgPCC.lasso.simulated.exp.v4(sample_size_set = my.sample_size_set,
                            gamma_set = rbind(c(0, 0, 0, 0, 0, 15, 10, 5, 5, 3, rep(0, p-10))),
							gamma_index = "3'",
                            rho = my.rho, 
                            p = 12, 
                            N = my.N, 
                            lambda_set = seq(0, 20, 1), 
                            time_units = "min", 
                            tol_0 = my.tol,
                            wd = mywd,
                            save.name="nonlead-new-algo")
}

# =================================================================================
# =================================================================================
# DIMENSION p = 30 (probably reduce number of replications?)
# =================================================================================
# =================================================================================

# =================================================================================
# =================================================================================
# DIMENSION p = 70 (probably reduce number of replications?)
# =================================================================================
# =================================================================================
















# =================================================================================
# sparsity 1, debug
# =================================================================================

if (casedebug == TRUE){
case1results <- RgPCC.lasso.simulated.exp.v4(sample_size_set = my.sample_size_set[1],
                          gamma_set = rbind(c(25, rep(0, p-1))),
			  gamma_index = "0",
                          rho = my.rho, 
                          p = p, 
                          N = 3, 
                          lambda_set = seq(0, 60, 2), 
                          time_units = "min", 
                          tol_0 = my.tol,
                          wd = mywd,
                          save.name="debug-new-algo")
}
