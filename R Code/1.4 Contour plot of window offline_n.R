library(plotly)
library(ggplot2)
options(scipen=999)  # turn off scientific notation like 1e+06
library(readr)
library(dplyr)
library(scales)
library(quantmod)
library(gridExtra)
library(grid)
library(webshot)

T=read.table("/Users/yangwen/Documents/Research/Programing/Data/Performance_MeanSqrt1_100_1d_3_n_Final23052019.csv",
#T=read.table("/Users/yangwen/Documents/Research/Programing/Data/Performance_Var2_n_23052019.csv",
             header = TRUE,
             sep = ",")
d = T$d
n = T$n
SAC = sqrt(T$AccuracyAC* T$SensitivityAC)
SMST = sqrt(T$AccuracyMST* T$SensitivityMST)
SR = sqrt(T$AccuracyR* T$SensitivityR)        

window = c(T$n[1], T$n[2], T$n[3], T$n[4], T$n[5], T$n[6], T$n[7],T$n[8])
dimension = (c(T$d[1],T$d[9], T$d[17], T$d[25], T$d[33],T$d[41],T$d[49]))
nR = length(window)
nC = length(dimension)
P_mean = matrix(c(SAC), nrow = nR, ncol = nC)

p1 = plot_ly(
  x = ~dimension, 
  y = ~window, 
  z = ~P_mean, 
  type = "contour",
  autocontour = F,
  contours = list(
    start = 0.1,
    end = 0.9,
    size = 0.05#, showlabels=TRUE 
  ),width = 1300, height = 400
) 
p1
p2 <- plot_ly(
  x = ~dimension, 
  y = ~window, 
  z = matrix(c(SMST), nrow = nR, ncol = nC), 
  type = "contour",#xaxis="Dimension (in natural log scale)", yaxis="P-mean",
  autocontour = F,
  contours = list(
    start = 0.1,
    end = 0.9,
    size = 0.05#, showlabels=TRUE
  ),width = 1300, height = 400
)
p2

p3 <- plot_ly(
  x = ~dimension, 
  y = ~window, 
  z = matrix(c(SR), nrow = nR, ncol = nC), 
  type = "contour",
  autocontour = F,
  contours = list(
    start = 0.1,
    end = 0.9,
    size = 0.05#, showlabels=TRUE
  ),width = 1300, height = 400
)
p3

p <- subplot(p1,p2,p3)
export(p, file = "/Users/yangwen/Desktop/Contour.png", selenium = NULL)

p4 = plot_ly(
  x = ~dimension, 
  y = ~window, 
  z = ~P_mean, 
  type = "contour",
  autocontour = F,
  contours = list(
    start = 0.1,
    end = 0.9,
    size = 0.05#, showlabels=TRUE
  ),width = 300, height = 600
) 
p4
export(p4, file = "/Users/yangwen/Desktop/Bar.png", selenium = NULL)
