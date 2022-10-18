## packages
library(MASS)
library(Rlab)
library(matlib)
library(xtable)
library(glmnet)
library(elasticnet)
library(caret)
library(grid)
library(gridExtra)
library(dplyr)
library(multcomp)
library(stargazer)
library(e1071)

# ========================================================================
# ========================================================================
# GRAPHICS
# ========================================================================
# ========================================================================


# =================================================================
# tune.visuals
# =================================================================
# the graphs for tuning process of RgPCC
# -----------------------------------------------------------------

tune.visuals <- function(AIC, BIC, MSE, pMSE, MSECV, myname = "", mycaption = ""){
  
  pdf(paste(myname, "-graphtuningAIC.pdf", sep = ""))
  myplotAIC <- ggplot(data = AIC, aes(x = lambda, y = AIC, group = sample_size)) +
    geom_line(aes(color = factor(sample_size))) +
    geom_point(aes(color = factor(sample_size))) +
    labs(title = "Parameter Tuning AIC", subtitle = mycaption, color = "Sample Size", y = "AIC")
  
  print(myplotAIC)
  dev.off()
  
  pdf(paste(myname, "-graphtuningBIC.pdf", sep = ""))
  myplotBIC <- ggplot(data = BIC, aes(x = lambda, y = BIC, group = sample_size)) +
    geom_line(aes(color = factor(sample_size))) +
    geom_point(aes(color = factor(sample_size))) +
    labs(title = "Parameter Tuning BIC", subtitle = mycaption, color = "Sample Size", y = "BIC")
  
  print(myplotBIC)
  dev.off()
  
  pdf(paste(myname, "-graphtuningMSE.pdf", sep = ""))
  myplotMSE <- ggplot(data = MSE, aes(x = lambda, y = MSE, group = sample_size)) +
    geom_line(aes(color = factor(sample_size))) +
    geom_point(aes(color = factor(sample_size))) +
    labs(title = "Parameter Tuning MSE", subtitle = mycaption, color = "Sample Size", y = "MSE")
  
  print(myplotMSE)
  dev.off()
  
  pdf(paste(myname, "-graphtuningpMSE.pdf", sep = ""))
  myplotpMSE <- ggplot(data = pMSE, aes(x = lambda, y = pMSE, group = sample_size)) +
    geom_line(aes(color = factor(sample_size))) +
    geom_point(aes(color = factor(sample_size))) +
    labs(title = "Parameter Tuning pMSE", subtitle = mycaption, color = "Sample Size", y = "pMSE")
  
  print(myplotpMSE)
  dev.off()
  
  pdf(paste(myname, "-graphtuningMSECV.pdf", sep = ""))
  myplotMSECV <- ggplot(data = MSECV, aes(x = lambda, y = MSECV, group = sample_size)) +
    geom_line(aes(color = factor(sample_size))) +
    geom_point(aes(color = factor(sample_size))) +
    labs(title = "Parameter Tuning MSECV", subtitle = mycaption, color = "Sample Size", y = "MSECV")
  
  print(myplotMSECV)
  dev.off()
  
}

latex.table <- function(results, save.name, mylabel = "", mycaption = "", myname = ""){
  capture.output(
    xtable(as.matrix(results),
           label = mylabel,
           caption = mycaption),
    file = paste(myname, ".tex", sep = ""),
    type = "output"
  )
}

metric.visuals <- function(results, save.name, mylabel = "", mycaption = "", myname = ""){
  
  #results <- results[order(results$method),]
  #results[,-1] <- round(as.numeric(results[,-1]), 4)
# =================================================================
# pdf table
# =================================================================
#
# -----------------------------------------------------------------
  g <- tableGrob(results)
  h <- grid::convertHeight(sum(g$heights), "mm", TRUE)
  w <- grid::convertHeight(sum(g$widths), "mm", TRUE)
  ggplot2::ggsave(paste(myname, ".pdf", sep = ""), g, height = h, width = w+3, units = "mm")






# =================================================================
# latex table
# =================================================================
# 
# -----------------------------------------------------------------

  capture.output(
    print(
      xtable(as.matrix(results),
             label = mylabel,
             caption = mycaption),
      file = paste(myname, ".tex", sep = ""),
      type = "latex"
    ),
    include.rownames = FALSE
  )
 
 
 
 
 
# =================================================================
# ratio of p's table
# =================================================================
# comments
#
#
#
# -----------------------------------------------------------------



}

