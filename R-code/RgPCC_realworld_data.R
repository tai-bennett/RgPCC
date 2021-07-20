# This script takes realworld data and compares the RgPCC method (tuned with AIC and BIC) and compares 
# the results to Logistic Regression, PC Logistic Regression and Elastic Net. The method "binary-classification-comparision.R"
# runs 5-fold CV to compare the errors of each method and performs a Tukey PW comparison test to see if the errors are 
# significantly different.

#laptop <- F

# Here, we source RgPCC, some functions and the binary classification method
#if (laptop == FALSE){
#  source("D:/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/RgPCC_lasso_experimental.R")
#  source("D:/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/RTG_functions.R")
#  source("D:/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/binary-classification-comparison.R")
#  setwd("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results")
#} else {
#  source("C:/Users/benne/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/RgPCC_lasso_experimental.R")
#  source("C:/Users/benne/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/RTG_functions.R")
#  source("C:/Users/benne/Dropbox/UA Documents/RTG_duncan_bennett/R_files_spring/binary-classification-comparison.R")
#  setwd("C:/Users/benne/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results")
#  }

  source('RgPCC_lasso_experimental_v2.R')
  source('RTG_functions.R')
  source('binary-classification-comparison.R')
  
library(tidyverse)
library(rlang)
library(readxl)
  

data4 = T # divorce data (YES but ties) : this results in a 3 way tie between AIC BIC and ENET
data9 = F # audit data (YES) : pretty good results
data12 = F # cryotherapy data (YES but ties) :
data14 = F # ecoli data (OK)

## ======================================================================================================
## 4 DIVORCE DATA
## ======================================================================================================


if (data4 == TRUE){
  divorce.file <- read.csv("data/divorce.csv", sep=";")
  
  set.seed(987)
  random.rows <- sample(nrow(divorce.file), 120)
  
  #divorce.all <- na.omit(divorce.file)[random.rows,]
  divorce.all <- na.omit(divorce.file)
  #divorce.all <- na.omit(divorce.file)
  rownames(divorce.all) <- NULL
  divorce.class <- divorce.all[55] # convert from b,g classes to 0,1 b = 0 and g = 1
  divorce.data <- divorce.all[-c(55)] # subset to get predictor data without id numbers and NAs
  divorce.data <- scale(divorce.data, center = TRUE, scale = FALSE) # center the data
  divorce.class.data <- cbind(divorce.data, divorce.class)
  divorce.class <- as.matrix(divorce.class)
  
  divorce.results <- bin.method.comp(
    X = divorce.data, 
    Y = divorce.class,
    lambda.set = seq(0, 5, 0.1),
    fold.size = 5,
    myseed = 120043,
    "divorce-5CV-errors", 
    wd = "data",
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "divorce-results.txt")
  print(divorce.results)
  sink(file = NULL)
}


## ======================================================================================================
## 9  AUDIT DATA
## ======================================================================================================

if (data9 == TRUE){
    audit.risk <- read.csv("data/audit_risk.csv")
	
audit.risk.all <- na.omit(audit.risk)

audit.risk.class <- audit.risk.all[27] # 0 = B and 1 = M
audit.risk.class <- as.matrix(audit.risk.class)

audit.risk.data <- audit.risk.all[-c(2, 5, 8, 12, 15, 18, 21, 23, 24, 25, 26, 27)] # subset to get predictor

audit.risk.data <- scale(audit.risk.data, center = TRUE, scale = FALSE) # center the data

audit.risk.results <- bin.method.comp(
  X = audit.risk.data, 
  Y = audit.risk.class,
  lambda.set = seq(0, 2, 0.1),
  fold.size = 5,
  myseed = 111,
  "audit.risk-5CV-errors", 
  wd = 'data',
  track = TRUE,
  print = TRUE
)


sink(file = "audit.risk-results.txt")
print(audit.risk.results)
sink(file = NULL)
}

## ======================================================================================================
## 12  CRYOTHERAPY DATA
## ======================================================================================================

if (data12 == TRUE){
  
    Cryotherapy <- read_excel("data/Cryotherapy.xlsx")

  Cryotherapy.all <- na.omit(Cryotherapy)
  
  Cryotherapy.class <- Cryotherapy.all[7]
  Cryotherapy.class <- as.matrix(Cryotherapy.class)
  
  Cryotherapy.data <- Cryotherapy.all[-c(7)] # subset to get predictor
  
  Cryotherapy.data <- scale(Cryotherapy.data, center = TRUE, scale = FALSE) # center the data
  
  Cryotherapy.results <- bin.method.comp(
    X = Cryotherapy.data, 
    Y = Cryotherapy.class,
    lambda.set = seq(0, 1, 0.01),
    fold.size = 5,
    myseed = 111,
    "Cryotherapy-5CV-errors", 
    wd = 'data',
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "Cryotherapy-results.txt")
  print(Cryotherapy.results)
  sink(file = NULL)
}




## ======================================================================================================
## 14  ECOLI DATA
## ======================================================================================================

if (data14 == TRUE){
  ecoli <- read.table("data/ecoli.data", row.names=1, quote="\"", comment.char="")
  ecoli.all <- na.omit(ecoli)
  
  ecoli.class <- I(ecoli.all[8] == "cp")
  ecoli.class <- as.matrix(ecoli.class)
  
  ecoli.data <- ecoli.all[-c(8)] # subset to get predictor
  
  ecoli.data <- scale(ecoli.data, center = TRUE, scale = FALSE) # center the data
  
  ecoli.results <- bin.method.comp(
    X = ecoli.data, 
    Y = ecoli.class,
    lambda.set = seq(0, 1, 0.1),
    fold.size = 5,
    myseed = 111,
    "ecoli-5CV-errors", 
    wd = "data",
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "ecoli-results.txt")
  print(ecoli.results)
  sink(file = NULL)
}


