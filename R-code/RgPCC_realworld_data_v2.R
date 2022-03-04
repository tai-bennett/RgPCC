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


library(tidyverse)
library(rlang)
library(readxl)
  
# Here we load some realworld data, do some cleaning and call the binary-classification-comparison method.
# YES : 1, 4, 9, 12, 114
# NO  : 2, 5, 8, 11, 15, 16

# which datasets to run?

data1 = F # the wisconsin breast cancer data (YES but lose to EN)
data2 = F # all fire data
data3 = F # coimbra
data4 = T # divorce data (YES but ties) : this results in a 3 way tie between AIC BIC and ENET +++++++++++
data5 = F # ionosphere
data6 = F # wpbc wisconsin prognostic breast cancer data
data7 = F # wdbc wisconsin diagnostic breast cancer data
data8 = F # hcc-survival data
data9 = F # audit data (YES) : pretty good results ++++++++++++++++++++
data10 = F # caesarian data : XXX 
data11 = F # ceramic data (YES) : pretty good but not great
data12 = F # cryotherapy data (YES but ties) : ++++++++++++++++++++++++
data13 = F
data14 = F # ecoli data (OK) +++++++++++++++++++++++++++++++++
data15 = F # enery data XXX
data16 = F # glass data
data17 = F # hcv egypt data


datanew = F
bigdata = F # for accumulation of data for 571B project


## ======================================================================================================
## 1 WISCONSIN BREAST CANCER DATA
## ======================================================================================================

# this data is subsetted but we may not have to do that

if (data1 == TRUE){
breast.cancer.wisconsin <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/breast-cancer-wisconsin/breast-cancer-wisconsin.txt", header=FALSE, na.strings="?")
names(breast.cancer.wisconsin) <- c("id", "clump.thickness", "unif.cell.size", "unif.cell.shape", "marginal.adhesion", "epithelial.cell.size", "bare.nuclei", "bland.chromatin", "normal.nucleoli", "mitoses", "class")

wisconsin.all <- na.omit(breast.cancer.wisconsin)

set.seed(120112)
random.rows <- sample(nrow(breast.cancer.wisconsin), 50)


wisconsin.all <- wisconsin.all[random.rows,]
rownames(wisconsin.all) <- NULL
wisconsin.class <- (wisconsin.all$class-2)/2 # convert from 2,4 classes to 0,1
wisconsin.data <- na.omit(wisconsin.all[-c(1,11)]) # subset to get predictor data without id numbers and NAs





wisconsin.data <- scale(wisconsin.data, center = TRUE, scale = FALSE) # center the data
wisconsin.class.data <- cbind(wisconsin.class, wisconsin.data)
wisconsin.class <- as.matrix(wisconsin.class)


wisconsin.results <- bin.method.comp(
  X = wisconsin.data, 
  Y = wisconsin.class,
  lambda.set = seq(0, 10, 0.1),
  fold.size = 5,
  myseed = 123,
  "wisconsin-5CV-errors", 
  wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
  track = TRUE
)


sink(file = "wisconsin-results.txt")
print(wisconsin.results)
sink(file = NULL)

dim(wisconsin.data)
}


## ======================================================================================================
## 2 FIRE DATA
## ======================================================================================================


if (data2 == TRUE){
Algerian.forest.fire.Bejaia <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/fire/Algerian-forest-fire-Bejaia.csv")

bejaia.fire.all <- na.omit(Algerian.forest.fire.Bejaia)
bejaia.fire.class <- I(bejaia.fire.all[14] == "fire") # convert not fire = 0 and fire = 1
bejaia.fire.data <- bejaia.fire.all[-c(1,2,3,14)] # subset to get predictor data without dates
bejaia.fire.data$days <- seq.int(nrow(bejaia.fire.data))-1 # these are days since june 1 2012

bejaia.fire.data <- scale(bejaia.fire.data, center = TRUE, scale = FALSE) # center the data
bejaia.fire.class.data <- cbind(bejaia.fire.data, bejaia.fire.class)
bejaia.fire.class <- as.matrix(bejaia.fire.class)

Algerian.forest.fire.Sidi.BelAbbes <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/fire/Algerian-forest-fire-Sidi-BelAbbes.csv")

Sidi.BelAbbes.fire.all <- na.omit(Algerian.forest.fire.Sidi.BelAbbes)
Sidi.BelAbbes.fire.class <- I(Sidi.BelAbbes.fire.all[14] == "fire") # convert not fire = 0 and fire = 1
Sidi.BelAbbes.fire.data <- Sidi.BelAbbes.fire.all[-c(1,2,3,14)] # subset to get predictor data without dates
Sidi.BelAbbes.fire.data$days <- seq.int(nrow(Sidi.BelAbbes.fire.data))-1 # these are days since june 1 2012

Sidi.BelAbbes.fire.data <- Sidi.BelAbbes.fire.data[-44,]
Sidi.BelAbbes.fire.class <- Sidi.BelAbbes.fire.class[-44,]

Sidi.BelAbbes.fire.data[, c(7,10)] <- sapply(Sidi.BelAbbes.fire.data[, c(7,10)], as.numeric)

Sidi.BelAbbes.fire.data <- scale(Sidi.BelAbbes.fire.data, center = TRUE, scale = FALSE) # center the data
Sidi.BelAbbes.fire.class.data <- cbind(Sidi.BelAbbes.fire.data, Sidi.BelAbbes.fire.class)
Sidi.BelAbbes.fire.class <- as.matrix(Sidi.BelAbbes.fire.class)

# combine these two data sets

all.fire.data <- rbind(bejaia.fire.data, Sidi.BelAbbes.fire.data)
all.fire.class <- rbind(bejaia.fire.class, Sidi.BelAbbes.fire.class)

all.fire.results <- bin.method.comp(
  X = all.fire.data, 
  Y = all.fire.class,
  lambda.set = seq(0, 15, 0.5),
  fold.size = 5,
  myseed = 173,
  "all.fire-5CV-errors", 
  wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
  track = TRUE
)


sink(file = "all-fire-results.txt")
print(all.fire.results)
sink(file = NULL)
}

## ======================================================================================================
## 3 COIMBRA DATA
## ======================================================================================================


if (data3 == TRUE){
coimbra <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/BreastCancerCoimbra/dataR2.csv")

coimbra.all <- na.omit(coimbra)

coimbra.class <- coimbra.all[10]
coimbra.data <- coimbra.all[-c(10)] # subset to get predictor

coimbra.data <- scale(coimbra.data, center = TRUE, scale = FALSE) # center the data
coimbra.class <- as.matrix(coimbra.class) - 1

coimbra.results <- bin.method.comp(
  X = coimbra.data, 
  Y = coimbra.class,
  lambda.set = seq(0, 40, 0.5),
  fold.size = 5,
  myseed = 111,
  "coimbra-5CV-errors", 
  wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
  track = TRUE,
  print = TRUE
)


sink(file = "coimbra-results.txt")
print(coimbra.results)
sink(file = NULL)
}


## ======================================================================================================
## 4 DIVORCE DATA
## ======================================================================================================


if (data4 == TRUE){
  divorce.file <- read.csv("/home/duncan/GitProjects/RgPCC/R-code/data/divorce.csv", sep=";")
  
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
    wd = "/home/duncan/GitProjects/RgPCC/R-code/data",
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "divorce-results.txt")
  print(divorce.results)
  sink(file = NULL)
}


## ======================================================================================================
## 5 IONOSPHERE DATA
## ======================================================================================================


if (data5 == TRUE){
ionosphere.file <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/ionosphere/ionosphere.txt", header=FALSE, na.strings="?")

set.seed(12019922)
random.rows <- sample(nrow(ionosphere.file), 55)

#ionosphere.all <- na.omit(ionosphere.file)[random.rows,]
ionosphere.all <- na.omit(ionosphere.file)
rownames(ionosphere.all) <- NULL
ionosphere.class <- I(ionosphere.all[35] == "g") # convert from b,g classes to 0,1 b = 0 and g = 1
ionosphere.data <- ionosphere.all[-c(2,35)] # subset to get predictor data without id numbers and NAs
ionosphere.data <- scale(ionosphere.data, center = TRUE, scale = FALSE) # center the data
ionosphere.class.data <- cbind(ionosphere.data, ionosphere.class)
ionosphere.class <- as.matrix(ionosphere.class)

ionosphere.results <- bin.method.comp(
  X = ionosphere.data, 
  Y = ionosphere.class,
  lambda.set = seq(0, 12, 0.5),
  fold.size = 5,
  myseed = 120043,
  "ionosphere-5CV-errors", 
  wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
  track = TRUE,
  print = TRUE
)

sink(file = "ionosphere-results.txt")
print(ionosphere.results)
sink(file = NULL)

}


## ======================================================================================================
## 6  WPBC DATA
## ======================================================================================================

if (data6 == TRUE){
  wpbc <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/breast-cancer-wisconsin/wpbc.txt", header=FALSE, na.strings = "?")
  
  wpbc.all <- na.omit(wpbc)
  
  wpbc.class <- I(wpbc.all[2] == "R") # 0 = N and 1 = R
  wpbc.class <- as.matrix(wpbc.class)
  
  wpbc.data <- wpbc.all[-c(1,2)] # subset to get predictor
  
  wpbc.data <- scale(wpbc.data, center = TRUE, scale = FALSE) # center the data
  
  wpbc.results <- bin.method.comp(
    X = wpbc.data, 
    Y = wpbc.class,
    lambda.set = seq(0, 40, 1),
    fold.size = 5,
    myseed = 111,
    "wpbc-5CV-errors", 
    wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "wpbc-results.txt")
  print(wpbc.results)
  sink(file = NULL)
  }

## ======================================================================================================
## 7  WDBC DATA
## ======================================================================================================

if (data7 == TRUE){
  wdbc <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/breast-cancer-wisconsin/wdbc.txt", header=FALSE, na.strings = "?")
  
  wdbc.all <- na.omit(wdbc)
  
  wdbc.class <- I(wdbc.all[2] == "M") # 0 = B and 1 = M
  wdbc.class <- as.matrix(wdbc.class)
  
  wdbc.data <- wdbc.all[-c(1,2)] # subset to get predictor
  
  wdbc.data <- scale(wdbc.data, center = TRUE, scale = FALSE) # center the data
  
  wdbc.results <- bin.method.comp(
    X = wdbc.data, 
    Y = wdbc.class,
    lambda.set = seq(0, 40, 1),
    fold.size = 5,
    myseed = 111,
    "wdbc-5CV-errors", 
    wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "wdbc-results.txt")
  print(wdbc.results)
  sink(file = NULL)
}


## ======================================================================================================
## 8  HCC-SURVIVAL DATA
## ======================================================================================================

## ======================================================================================================
## 9  AUDIT DATA
## ======================================================================================================

if (data9 == TRUE){
  if (laptop == FALSE) {
    audit.risk <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/audit_data/audit_risk.csv")
    mywd <- "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results"
  } else {
    audit.risk <- read.csv("C:/Users/benne/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/audit_data/audit_risk.csv")
    mywd <- "C:/Users/benne/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results"
  }
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
  wd = mywd,
  track = TRUE,
  print = TRUE
)


sink(file = "audit.risk-results.txt")
print(audit.risk.results)
sink(file = NULL)
}

## ======================================================================================================
## 10  CAESARIAN SETION DATA
## ======================================================================================================

if (data10 == TRUE){
  caesarian.section <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/caesarian-section.csv")
  caesarian.section.all <- na.omit(caesarian.section)
  
  caesarian.section.class <- caesarian.section.all[7] # 0 = B and 1 = M
  caesarian.section.class <- as.matrix(caesarian.section.class)
  
  caesarian.section.data <- caesarian.section.all[-c(1,7)] # subset to get predictor
  
  caesarian.section.data <- scale(caesarian.section.data, center = TRUE, scale = FALSE) # center the data
  
  caesarian.section.results <- bin.method.comp(
    X = caesarian.section.data, 
    Y = caesarian.section.class,
    lambda.set = seq(0, 3, 0.05),
    fold.size = 5,
    myseed = 1311,
    "caesarian.section-5CV-errors", 
    wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "caesarian.section-results.txt")
  print(caesarian.section.results)
  sink(file = NULL)
}

## ======================================================================================================
## 11  CERAMIC DATA
## ======================================================================================================

if (data11 == TRUE){
  Chemical.Composion.Ceramic <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/Chemical-Composion-Ceramic.csv", header=TRUE)
  Chemical.Composion.Ceramic$class <- I(startsWith(Chemical.Composion.Ceramic$Ceramic.Name, "DY"))
  Chemical.Composion.Ceramic$partnew <- I(Chemical.Composion.Ceramic$Part == "Body")
  
  Chemical.Composion.Ceramic.all <- na.omit(Chemical.Composion.Ceramic)
  
  Chemical.Composion.Ceramic.class <- Chemical.Composion.Ceramic.all[20] # 0 = B and 1 = M
  Chemical.Composion.Ceramic.class <- as.matrix(Chemical.Composion.Ceramic.class)
  
  Chemical.Composion.Ceramic.data <- Chemical.Composion.Ceramic.all[-c(1, 2, 20)] # subset to get predictor
  
  Chemical.Composion.Ceramic.data <- scale(Chemical.Composion.Ceramic.data, center = TRUE, scale = FALSE) # center the data
  
  Chemical.Composion.Ceramic.results <- bin.method.comp(
    X = Chemical.Composion.Ceramic.data, 
    Y = Chemical.Composion.Ceramic.class,
    lambda.set = seq(0, 2, 0.01),
    fold.size = 5,
    myseed = 111,
    "Chemical.Composion.Ceramic-5CV-errors", 
    wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "Chemical.Composion.Ceramic-results.txt")
  print(Chemical.Composion.Ceramic.results)
  sink(file = NULL)
}

## ======================================================================================================
## 12  CRYOTHERAPY DATA
## ======================================================================================================

if (data12 == TRUE){
  
  if (laptop == FALSE) {
    Cryotherapy <- read_excel("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/Cryotherapy.xlsx")
    mywd <- "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results"
  } else {
    Cryotherapy <- read_excel("C:/Users/benne/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/Cryotherapy.xlsx")
    mywd <- "C:/Users/benne/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results"
  }

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
    wd = mywd,
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "Cryotherapy-results.txt")
  print(Cryotherapy.results)
  sink(file = NULL)
}

## ======================================================================================================
## 13 MESSIDOR (DIABETES) DATA
## ======================================================================================================

if (data13 == TRUE){
  messidor.features <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/messidor_features.csv")
  messidor.features.all <- na.omit(messidor.features)
  
  messidor.features.class <- messidor.features.all[19] # 0 = B and 1 = M
  messidor.features.class <- as.matrix(messidor.features.class)
  
  messidor.features.data <- messidor.features.all[-c(1, 2, 3, 20, 21)] # subset to get predictor
  
  messidor.features.data <- scale(messidor.features.data, center = TRUE, scale = FALSE) # center the data
  
  messidor.features.results <- bin.method.comp(
    X = messidor.features.data, 
    Y = messidor.features.class,
    lambda.set = seq(0, 10, 1),
    fold.size = 5,
    myseed = 111,
    "messidor.features-5CV-errors", 
    wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "messidor.features-results.txt")
  print(messidor.features.results)
  sink(file = NULL)
}




## ======================================================================================================
## 14  ECOLI DATA
## ======================================================================================================

if (data14 == TRUE){
  ecoli <- read.table("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/ecoli.data", row.names=1, quote="\"", comment.char="")
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
    wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "ecoli-results.txt")
  print(ecoli.results)
  sink(file = NULL)
}


## ======================================================================================================
## 15  ENERGY DATA
## ======================================================================================================

if (data15 == TRUE){
  energy <- read_excel("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/energy.xlsx")
  energy.all <- na.omit(energy)
  
  energy.class <- I(energy.all[9] > 15) # 0 = B and 1 = M
  energy.class <- as.matrix(energy.class)
  
  energy.class.cool <- I(energy.all[10] > 15) # 0 = B and 1 = M
  energy.class.cool <- as.matrix(energy.class)
  
  energy.data <- energy.all[-c(9, 10)] # subset to get predictor
  
  energy.data <- scale(energy.data, center = TRUE, scale = FALSE) # center the data
  
  energy.results <- bin.method.comp(
    X = energy.data, 
    Y = energy.class,
    lambda.set = seq(0, 0.1, 0.001),
    fold.size = 5,
    myseed = 118751,
    "energy-5CV-errors", 
    wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "energy-results.txt")
  print(energy.results)
  sink(file = NULL)
}

## ======================================================================================================
## 16  GLASS DATA
## ======================================================================================================

if (data16 == TRUE){
  glass <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/glass.data", header=FALSE, row.names=1)
  glass.all <- na.omit(glass)
  
  glass.class <- I(glass.all[10] == 1) # 0 = B and 1 = M
  glass.class <- as.matrix(glass.class)
  
  glass.data <- glass.all[-c(10)] # subset to get predictor
  
  glass.data <- scale(glass.data, center = TRUE, scale = FALSE) # center the data
  
  glass.results <- bin.method.comp(
    X = glass.data, 
    Y = glass.class,
    lambda.set = seq(0, 5, 0.1),
    fold.size = 5,
    myseed = 111,
    "glass-5CV-errors", 
    wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "glass-results.txt")
  print(glass.results)
  sink(file = NULL)
}

## ======================================================================================================
## 17 HCV EGYPT
## ======================================================================================================

if (data17 == TRUE){
  HCV.Egy.Data <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/HCV-Egy-Data.csv", header=TRUE)
  glass.all <- na.omit(glass)
  
  glass.class <- I(glass.all[10] == 1) # 0 = B and 1 = M
  glass.class <- as.matrix(glass.class)
  
  glass.data <- glass.all[-c(10)] # subset to get predictor
  
  glass.data <- scale(glass.data, center = TRUE, scale = FALSE) # center the data
  
  glass.results <- bin.method.comp(
    X = glass.data, 
    Y = glass.class,
    lambda.set = seq(0, 5, 0.1),
    fold.size = 5,
    myseed = 111,
    "glass-5CV-errors", 
    wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "glass-results.txt")
  print(glass.results)
  sink(file = NULL)
}



## ======================================================================================================
## ======================================================================================================
## error data for 571B analysis
## ======================================================================================================
## ======================================================================================================


if (bigdata == TRUE) {
error.dataset.credit <- credit.results$error.results
error.dataset.credit$data <- 1

error.dataset.ionosphere <- ionosphere.results$error.results
error.dataset.ionosphere$data <- 2

error.dataset.wpbc <- wpbc.results$error.results
error.dataset.wpbc$data <- 3


error.dataset.coimbra <- coimbra.results$error.results
error.dataset.coimbra$data <- 4
}

#all.dataset <- rbind(error.dataset.credit, error.dataset.ionosphere, error.dataset.wpbc, error.dataset.coimbra)
#colnames(all.dataset) <- c("cverror", "modeltype", "data")
#all.dataset <- all.dataset[order(all.dataset$data),]
#all.dataset <- all.dataset[,c(1,3,2)]
#write.csv(all.dataset, "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results/alldataset.csv", row.names = T)



#error.data.fit <- lm(cverror ~ factor(modeltype) + factor(data) + factor(modeltype)*factor(cverror), data = all.dataset)
#eerror.resid <- residuals(error.data.fit)
#error.fitted <- fitted(error.data.fit)

#bc.res <- boxcox(error.data.fit, lambda = seq(-2, 2, 0.1), inter = F)


#sasLM::GLM(cverror ~ factor(modeltype) + factor(data) + factor(modeltype)*factor(cverror), all.dataset)




## ======================================================================================================
## ??  NEW DATA
## ======================================================================================================

if (datanew == TRUE){
  new <- read.csv("D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/audit_data/audit_risk.csv")
  new.all <- na.omit(new)
  
  new.class <- new.all[1] # 0 = B and 1 = M
  new.class <- as.matrix(new.class)
  
  new.data <- new.all[-c(1)] # subset to get predictor
  
  new.data <- scale(new.data, center = TRUE, scale = FALSE) # center the data
  
  new.results <- bin.method.comp(
    X = new.data, 
    Y = new.class,
    lambda.set = seq(0, 10, 1),
    fold.size = 5,
    myseed = 111,
    "new-5CV-errors", 
    wd = "D:/Dropbox/UA Documents/RTG_duncan_bennett/RealworldData/error-results",
    track = TRUE,
    print = TRUE
  )
  
  
  sink(file = "new-results.txt")
  print(new.results)
  sink(file = NULL)
}

