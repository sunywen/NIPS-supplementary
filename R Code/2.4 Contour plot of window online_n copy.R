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

T=read.table("/Users/yangwen/Documents/Research/Programing/Data/Performance_2.2.1_ChgMeanB100_3rt_n_FPR_crop20052019.csv",
             header = TRUE,
             sep = ",")
H=read.table("/Users/yangwen/Documents/Research/Programing/Data/Performance_2.3.1_VarChgB100_2s_n_FPR_crop20052019.csv",
             header = TRUE,
             sep = ",")
d = T$d
n = T$n
SAC = sqrt(T$AccuracyAC* T$SensitivityAC)
SMST = sqrt(T$AccuracyMST* T$SensitivityMST)

VSAC = sqrt(H$AccuracyAC* H$SensitivityAC)
VSMST = sqrt(H$AccuracyMST* H$SensitivityMST)


window = c(T$n[1], T$n[2], T$n[3], T$n[4], T$n[5], T$n[6], T$n[7],T$n[8])
dimension = log(c(T$d[1],T$d[9], T$d[17], T$d[25], T$d[33],T$d[41]))
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
  ),width = 900, height = 400
) 
p1
p2 <- plot_ly(
  x = ~dimension, 
  y = ~window, 
  z = matrix(c(SMST), nrow = nR, ncol = nC), 
  type = "contour",
  autocontour = F,
  contours = list(
    start = 0.1,
    end = 0.9,
    size = 0.05#, showlabels=TRUE
  ),width = 900, height = 400
)
p2

p3 = plot_ly(
  x = ~dimension, 
  y = ~window, 
  z = matrix(c(VSAC), nrow = nR, ncol = nC), 
  type = "contour",
  autocontour = F,
  contours = list(
    start = 0.1,
    end = 0.9,
    size = 0.05#, showlabels=TRUE 
  ),width = 900, height = 400
) 

p4 <- plot_ly(
  x = ~dimension, 
  y = ~window, 
  z = matrix(c(VSMST), nrow = nR, ncol = nC), 
  type = "contour",
  autocontour = F,
  contours = list(
    start = 0.1,
    end = 0.9,
    size = 0.05#, showlabels=TRUE
  ),width = 900, height = 400
)
p2

p <- subplot(p3,p4)



export(p, file = "/Users/yangwen/Desktop/Contour.png", selenium = NULL)

p5 = plot_ly(
  x = ~dimension, 
  y = ~window, 
  z = ~P_mean, 
  type = "contour",
  autocontour = F,
  contours = list(
    start = 0.1,
    end = 0.9,
    size = 0.05, showlabels=TRUE
  ),width = 300, height = 600
) 
p5
export(p5, file = "/Users/yangwen/Desktop/Bar.png", selenium = NULL)
