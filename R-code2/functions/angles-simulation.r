library(MASS)
library(Rlab)
library(matlib)
library(ggplot2)
library(stargazer)
library(grid)
library(gridExtra)
library(xtable)

angle_simulation <- function(
    gamma_set,
    gamma_index,
    rho,
    p,
    sample_size_factor = 5,
    lambda_set,
    time_units = "min",
    tol_0 = 0.1,
    wd
) {
    # load in functions from the fuctino folder
    # assumes function is called from scripts folder
    current_dir  <- getwd()
    setwd("../functions")

    functions_list <-
    list.files(c("."), pattern = "*.r$", full.names = TRUE, ignore.case = TRUE)

    functions_list <-
    functions_list[!functions_list %in% c("./simulation_function.r")]

    sapply(functions_list, source, .GlobalEnv)
    setwd("../")
    setwd(current_dir)
}
