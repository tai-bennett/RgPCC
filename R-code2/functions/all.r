current_dir  <- getwd()

setwd("../functions")

functions_list <-
list.files(c("."), pattern = "*.r$", full.names = TRUE, ignore.case = TRUE)

functions_list <- functions_list[!functions_list %in% c("./all.r")]

sapply(functions_list, source, .GlobalEnv)

setwd(current_dir)
