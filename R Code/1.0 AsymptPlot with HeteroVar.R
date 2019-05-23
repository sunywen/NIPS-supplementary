options(scipen=999)  # turn off scientific notation like 1e+06
library(ggplot2)
library(readr)
library(dplyr)
library(scales)
library(quantmod)
library(gridExtra)
library(grid)

CV_table = read.table("/Users/yangwen/Desktop/Research/Programing/Data/Asymptotics_Heterogeneous.csv",
                      header = TRUE,
                      sep = ",")
CV_table$n=CV_table$n/2

g1 = ggplot(subset(CV_table,d %in% c(1)), aes(x = n, y =  CU)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g2 = ggplot(subset(CV_table,d %in% c(10)), aes(x = n, y =  CU)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g3 = ggplot(subset(CV_table,d %in% c(50)), aes(x = n, y =  CU)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g4 = ggplot(subset(CV_table,d %in% c(100)), aes(x = n, y =  CU)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g5 = ggplot(subset(CV_table,d %in% c(500)), aes(x = n, y =  CU)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g6 = ggplot(subset(CV_table,d %in% c(1)), aes(x = n, y =  CU2)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g7 = ggplot(subset(CV_table,d %in% c(10)), aes(x = n, y =  CU2)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g8 = ggplot(subset(CV_table,d %in% c(50)), aes(x = n, y =  CU2)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g9 = ggplot(subset(CV_table,d %in% c(100)), aes(x = n, y =  CU2)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g10 = ggplot(subset(CV_table,d %in% c(500)), aes(x = n, y =  CU2)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g11 = ggplot(subset(CV_table,d %in% c(1)), aes(x = n, y =  ACU)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g12 = ggplot(subset(CV_table,d %in% c(10)), aes(x = n, y =  ACU)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g13 = ggplot(subset(CV_table,d %in% c(50)), aes(x = n, y =  ACU)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g14 = ggplot(subset(CV_table,d %in% c(100)), aes(x = n, y =  ACU)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g15 = ggplot(subset(CV_table,d %in% c(500)), aes(x = n, y =  ACU)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g16 = ggplot(subset(CV_table,d %in% c(1)), aes(x = n, y =  ACU2)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g17 = ggplot(subset(CV_table,d %in% c(10)), aes(x = n, y =  ACU2)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g18 = ggplot(subset(CV_table,d %in% c(50)), aes(x = n, y =  ACU2)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g19 = ggplot(subset(CV_table,d %in% c(100)), aes(x = n, y =  ACU2)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))

g20 = ggplot(subset(CV_table,d %in% c(500)), aes(x = n, y =  ACU2)) + geom_point(aes(col = d), size=3) + labs(y = expression(rho[alpha]))



pdf("/Users/yangwen/Desktop/Research/Analysis/Asymptotics of the Tstat of mean - MST, Heterogeneous variance.pdf")

grid.arrange(
  g1,
  g2,
  g3,
  g4,
  g5,
  nrow = 3, ncol = 2,
  top = "Asymptotics of the Tstat of mean - MST, Heterogeneous variance", 
  bottom = textGrob(
    "2.1",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  ))
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Asymptotics of the Tstat of variance - MST, Heterogeneous variance.pdf")

grid.arrange(
  g6,
  g7,
  g8,
  g9,
  g10,
  nrow = 3, ncol = 2,
  top = "Asymptotics of the Tstat of variance - MST, Heterogeneous variance", 
  bottom = textGrob(
    "2.2",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  ))
dev.off()


pdf("/Users/yangwen/Desktop/Research/Analysis/Asymptotics of the Tstat of mean - AC, Heterogeneous variance.pdf")

grid.arrange(
  g11,
  g12,
  g13,
  g14,
  g15,
  nrow = 3, ncol = 2,
  top = "Asymptotics of the Tstat of mean - AC, Heterogeneous variance", 
  bottom = textGrob(
    "2.3",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  ))
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Asymptotics of the Tstat of variance - AC, Heterogeneous variance.pdf")

grid.arrange(
  g16,
  g17,
  g18,
  g19,
  g20,
  nrow = 3, ncol = 2,
  top = "Asymptotics of the Tstat of variance - AC, Heterogeneous variance", 
  bottom = textGrob(
    "2.4",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  ))
dev.off()
