# testing the new RgPCC function for simulated data

#source("D:/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/RgPCC_lasso_simulated_experiment_function.R")
#source("D:/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/RgPCC_lasso_simulated_experiment_function_v2.R")

#source('RgPCC_lasso_simulated_experiment_function.R')
source('RgPCC_lasso_simulated_experiment_function_v2.R')



casedebug = F

#lead cases
case1 = T #
case2 = T #
case3 = F #
case4 = T #

#nonlead cases
case5 = F #
case6 = T #
case7 = T #
case8 = T #

p <- 80 # number of predictor variables
my.sample_size_set <- c(200,300) # sample size of each sample

# set of true gammas
my.gamma_set <- rbind(
  c(25, rep(0, p-1)),
  c(20, 10, 10, rep(0, p-3)),
  c(15, 10, 5, 5, 3, rep(0, p-5)),
  c(20, 0, 0.4, 0.8, 0.8, 1.2, 0.4, rep(0, p-7)),
  c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, rep(0, p-12)),
  c(0, 0, 0.4, 0.8, 0.8, 1.2, 0.4, 0, 0, 0, 0, 20, rep(0, p-12))
)

# covariance matrix generation
my.rho <- 0.8


# ------------------------------------- parameters to vary over each test
# tolerance of Newton-Ralphson method
m.ytol_0 <- 0.1

# tuning parameters and sample sizes
my.lambda_set <- c(50)

my.sample_size_factor <- 5 # how much larger the testing set is than the training set

my.N <- 100 # N is number of samples


# =================================================================================
# =================================================================================
# DIMENSION p = 12
# =================================================================================
# =================================================================================

# =================================================================================
# sparsity 1, N = 100, p = 12 (1)
# =================================================================================

if (case1 == TRUE){
RgPCC.lasso.simulated.exp.v2(sample_size_set = c(200),
                          gamma_set = rbind(c(25, rep(0, p-1))),
                          rho = 0.8, 
                          p = p, 
                          N = 100, 
                          lambda_set = seq(0, 60, 2), 
                          time_units = "sec", 
                          tol_0 = 0.1,
                          wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/Simulated_Experiments_v2",
                          save.name="sparsity1-100")
}

# =================================================================================
# sparsity 1, N = 200, p = 12 (2)
# =================================================================================

if (case2 == TRUE){
RgPCC.lasso.simulated.exp.v2(sample_size_set = c(300),
                          gamma_set = rbind(c(25, rep(0, p-1))),
                          rho = 0.8, 
                          p = p, 
                          N = 100, 
                          lambda_set = seq(0, 150, 5), 
                          time_units = "sec", 
                          tol_0 = 0.1,
                          wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/Simulated_Experiments_v2",
                          save.name="sparsity1-200")
}

# =================================================================================
# sparsity 3, N = {100, 200}, p = 12
# =================================================================================

if (case3 == TRUE){
RgPCC.lasso.simulated.exp.v2(sample_size_set = c(200,300),
                          gamma_set = rbind(c(25, rep(0, p-1))),
                          rho = 0.8, 
                          p = p, 
                          N = 100, 
                          lambda_set = seq(0, 40, 2), 
                          time_units = "sec", 
                          tol_0 = 0.1,
                          wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/Simulated_Experiments_v2",
                          save.name="sparsity3")
}

# =================================================================================
# sparsity 5, N = 100,200, p = 12
# =================================================================================

if (case4 == TRUE){
RgPCC.lasso.simulated.exp.v2(sample_size_set = c(200,300),
                          gamma_set = rbind(c(15, 10, 5, 5, 3, rep(0, p-5))),
                          rho = 0.8, 
                          p = p, 
                          N = 100, 
                          lambda_set = seq(0, 60, 2), 
                          time_units = "sec", 
                          tol_0 = 0.1,
                          wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/Simulated_Experiments_v2",
                          save.name="sparsity5")
}

# =================================================================================
# =================================================================================
# Nonleading sparsity
# =================================================================================
# =================================================================================

# =================================================================================
# sparsity 1, N = 100, p = 12 (1)
# =================================================================================

if (case5 == TRUE){
  RgPCC.lasso.simulated.exp.v2(sample_size_set = c(200),
                            gamma_set = rbind(c(0, 25, rep(0, p-2))),
                            rho = 0.8, 
                            p = p, 
                            N = 100, 
                            lambda_set = seq(0, 60, 2), 
                            time_units = "sec", 
                            tol_0 = 0.1,
                            wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/Simulated_Experiments_v2",
                            save.name="sparsity1-100-nonlead")
}


# =================================================================================
# sparsity 1, N = 100, p = 12
# =================================================================================

if (case6 == TRUE){
  RgPCC.lasso.simulated.exp.v2(sample_size_set = c(200,300),
                            gamma_set = rbind(c(0, 0, 0, 0, 0, 25, rep(0, p-6))),
                            rho = 0.8, 
                            p = p, 
                            N = 100, 
                            lambda_set = seq(0, 100, 10), 
                            time_units = "sec", 
                            tol_0 = 0.1,
                            wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/Simulated_Experiments_v2",
                            save.name="sparsity1-nonlead")
}


# =================================================================================
# sparsity 3, N = 100, p = 12 (1)
# c(20, 10, 10, rep(0, p-3)),
# c(15, 10, 5, 5, 3, rep(0, p-5))
# =================================================================================

if (case7 == TRUE){
  RgPCC.lasso.simulated.exp.v2(sample_size_set = c(200,300),
                            gamma_set = rbind(c(0, 0, 0, 0, 0, 20, 10, 10, rep(0, p-8))),
                            rho = 0.8, 
                            p = p, 
                            N = 100, 
                            lambda_set = seq(0, 20, 1), 
                            time_units = "sec", 
                            tol_0 = 0.1,
                            wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/Simulated_Experiments_v2",
                            save.name="sparsity3-nonlead")
}


# =================================================================================
# sparsity 5, N = 100, p = 12 (1)
# =================================================================================

if (case8 == TRUE){
  RgPCC.lasso.simulated.exp.v2(sample_size_set = c(200,300),
                            gamma_set = rbind(c(0, 0, 0, 0, 0, 15, 10, 5, 5, 3, rep(0, p-10))),
                            rho = 0.8, 
                            p = 12, 
                            N = 100, 
                            lambda_set = seq(0, 20, 1), 
                            time_units = "sec", 
                            tol_0 = 0.1,
                            wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/Simulated_Experiments_v2",
                            save.name="sparsity5-nonlead")
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
  RgPCC.lasso.simulated.exp.v2(sample_size_set = c(200,300),
                            gamma_set = rbind(c(25, rep(0, p-1))),
                            rho = 0.8, 
                            p = 12, 
                            N = 5, 
                            lambda_set = seq(0, 4, 2), 
                            time_units = "sec", 
                            tol_0 = 0.1,
                            wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/Simulated_Experiments",
                            save.name="debug")
}

