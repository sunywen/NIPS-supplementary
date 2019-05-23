options(scipen=999)  # turn off scientific notation like 1e+06
library(ggplot2)
library(readr)
library(dplyr)
library(scales)
library(quantmod)
library(gridExtra)
library(grid)

Power_table = read.table("/Users/yangwen/Desktop/Research/Programing/Data/Performance_Var2_100_2.csv",
                         header = TRUE,
                         sep = ",")
Power_table$n = Power_table$n/2

g1 = ggplot(subset(Power_table,d %in% c(1)), aes(n)) +   geom_point(aes(y = AccuracyMST, colour = "MST, d=1")) + 
  geom_point(aes(y = AccuracyAC, colour = "AC, d=1")) +
  geom_point(aes(y = AccuracyR, colour = "Ref.Chen, d=1"))

g2 = ggplot(subset(Power_table,d %in% c(10)), aes(n)) +   geom_point(aes(y = AccuracyMST, colour = "MST, d=10")) + 
  geom_point(aes(y = AccuracyAC, colour = "AC, d=10")) +
  geom_point(aes(y = AccuracyR, colour = "Ref.Chen, d=10"))

g3 = ggplot(subset(Power_table,d %in% c(50)), aes(n)) +   geom_point(aes(y = AccuracyMST, colour = "MST, d=50")) + 
  geom_point(aes(y = AccuracyAC, colour = "AC, d=50")) +
  geom_point(aes(y = AccuracyR, colour = "Ref.Chen, d=50"))

g4 = ggplot(subset(Power_table,d %in% c(100)), aes(n)) +   geom_point(aes(y = AccuracyMST, colour = "MST, d=100")) + 
  geom_point(aes(y = AccuracyAC, colour = "AC, d=100")) +
  geom_point(aes(y = AccuracyR, colour = "Ref.Chen, d=100"))

g5 = ggplot(subset(Power_table,d %in% c(500)), aes(n)) +   geom_point(aes(y = AccuracyMST, colour = "MST, d=500")) + 
  geom_point(aes(y = AccuracyAC, colour = "AC, d=500")) +
  geom_point(aes(y = AccuracyR, colour = "Ref.Chen, d=500"))

pdf("/Users/yangwen/Desktop/Research/Analysis/Test Accuracy for change in Variance.pdf")

grid.arrange(
  g1,
  g2,
  g3,
  g4,
  g5,
  nrow = 3, ncol = 2,
  top = "Test Accuracy- for change in Variance", 
  bottom = textGrob(
    "4.1",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  ))
dev.off()


g6 = ggplot(subset(Power_table,d %in% c(1)), aes(n)) +   geom_point(aes(y = SensitivityMST, colour = "MST, d=1")) + 
  geom_point(aes(y = SensitivityAC, colour = "AC, d=1")) +
  geom_point(aes(y = SensitivityR, colour = "Ref.Chen, d=1"))

g7 = ggplot(subset(Power_table,d %in% c(10)), aes(n)) +   geom_point(aes(y = SensitivityMST, colour = "MST, d=10")) + 
  geom_point(aes(y = SensitivityAC, colour = "AC, d=10")) +
  geom_point(aes(y = SensitivityR, colour = "Ref.Chen, d=10"))

g8 = ggplot(subset(Power_table,d %in% c(50)), aes(n)) +   geom_point(aes(y = SensitivityMST, colour = "MST, d=50")) + 
  geom_point(aes(y = SensitivityAC, colour = "AC, d=50")) +
  geom_point(aes(y = SensitivityR, colour = "Ref.Chen, d=50"))

g9 = ggplot(subset(Power_table,d %in% c(100)), aes(n)) +   geom_point(aes(y = SensitivityMST, colour = "MST, d=100")) + 
  geom_point(aes(y = SensitivityAC, colour = "AC, d=100")) +
  geom_point(aes(y = SensitivityR, colour = "Ref.Chen, d=100"))

g10 = ggplot(subset(Power_table,d %in% c(500)), aes(n)) +   geom_point(aes(y = SensitivityMST, colour = "MST, d=500")) + 
  geom_point(aes(y = SensitivityAC, colour = "AC, d=500")) +
  geom_point(aes(y = SensitivityR, colour = "Ref.Chen, d=500"))

pdf("/Users/yangwen/Desktop/Research/Analysis/Test Sensitivity for change in Variance.pdf")

grid.arrange(
  g6,
  g7,
  g8,
  g9,
  g10,
  nrow = 3, ncol = 2,
  top = "Test Sensitivity - for change in Variance", 
  bottom = textGrob(
    "4.2",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  ))
dev.off()