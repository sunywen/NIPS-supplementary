options(scipen=999)  # turn off scientific notation like 1e+06
library(ggplot2)
library(readr)
library(dplyr)
library(scales)
library(quantmod)
library(gridExtra)
library(grid)

CV_table = read.table("/Users/yangwen/Desktop/Research/Asymptotics_Homogeneous.csv",
                      header = TRUE,
                      sep = ",")
g1 = ggplot(subset(CV_table,d %in% c(1)), aes(x = n, y =  CU)) + geom_point(aes(col = d), size=3) 

g2 = ggplot(subset(CV_table,d %in% c(10)), aes(x = n, y =  CU)) + geom_point(aes(col = d), size=3) 

g3 = ggplot(subset(CV_table,d %in% c(50)), aes(x = n, y =  CU)) + geom_point(aes(col = d), size=3) 

g4 = ggplot(subset(CV_table,d %in% c(100)), aes(x = n, y =  CU)) + geom_point(aes(col = d), size=3) 

g5 = ggplot(subset(CV_table,d %in% c(500)), aes(x = n, y =  CU)) + geom_point(aes(col = d), size=3) 

g6 = ggplot(subset(CV_table,d %in% c(1)), aes(x = n, y =  CU2)) + geom_point(aes(col = d), size=3) 

g7 = ggplot(subset(CV_table,d %in% c(10)), aes(x = n, y =  CU2)) + geom_point(aes(col = d), size=3) 

g8 = ggplot(subset(CV_table,d %in% c(50)), aes(x = n, y =  CU2)) + geom_point(aes(col = d), size=3) 

g9 = ggplot(subset(CV_table,d %in% c(100)), aes(x = n, y =  CU2)) + geom_point(aes(col = d), size=3) 

g10 = ggplot(subset(CV_table,d %in% c(500)), aes(x = n, y =  CU2)) + geom_point(aes(col = d), size=3) 

g11 = ggplot(subset(CV_table,d %in% c(1)), aes(x = n, y =  ACU)) + geom_point(aes(col = d), size=3) 

g12 = ggplot(subset(CV_table,d %in% c(10)), aes(x = n, y =  ACU)) + geom_point(aes(col = d), size=3) 

g13 = ggplot(subset(CV_table,d %in% c(50)), aes(x = n, y =  ACU)) + geom_point(aes(col = d), size=3) 

g14 = ggplot(subset(CV_table,d %in% c(100)), aes(x = n, y =  ACU)) + geom_point(aes(col = d), size=3) 

g15 = ggplot(subset(CV_table,d %in% c(500)), aes(x = n, y =  ACU)) + geom_point(aes(col = d), size=3) 

g16 = ggplot(subset(CV_table,d %in% c(1)), aes(x = n, y =  ACU2)) + geom_point(aes(col = d), size=3) 

g17 = ggplot(subset(CV_table,d %in% c(10)), aes(x = n, y =  ACU2)) + geom_point(aes(col = d), size=3) 

g18 = ggplot(subset(CV_table,d %in% c(50)), aes(x = n, y =  ACU2)) + geom_point(aes(col = d), size=3) 

g19 = ggplot(subset(CV_table,d %in% c(100)), aes(x = n, y =  ACU2)) + geom_point(aes(col = d), size=3) 

g20 = ggplot(subset(CV_table,d %in% c(500)), aes(x = n, y =  ACU2)) + geom_point(aes(col = d), size=3) 


pdf("/Users/yangwen/Desktop/Research/Asymptotics of the Tstat of mean - MST, Homogeneous variance.pdf")

grid.arrange(
  g1,
  g2,
  g3,
  g4,
  g5,
  nrow = 3, ncol = 2,
  top = "Asymptotics of the Tstat of mean - MST, Homogeneous variance", 
  bottom = textGrob(
    "this footnote is right-justified",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  ))
dev.off()

pdf("/Users/yangwen/Desktop/Research/Asymptotics of the Tstat of variance - MST, Homogeneous variance.pdf")

grid.arrange(
  g6,
  g7,
  g8,
  g9,
  g10,
  nrow = 3, ncol = 2,
  top = "Asymptotics of the Tstat of variance - MST, Homogeneous variance", 
  bottom = textGrob(
    "this footnote is right-justified",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  ))
dev.off()


pdf("/Users/yangwen/Desktop/Research/Asymptotics of the Tstat of mean - AC, Homogeneous variance.pdf")

grid.arrange(
  g11,
  g12,
  g13,
  g14,
  g15,
  nrow = 3, ncol = 2,
  top = "Asymptotics of the Tstat of mean - AC, Homogeneous variance", 
  bottom = textGrob(
    "this footnote is right-justified",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  ))
dev.off()

pdf("/Users/yangwen/Desktop/Research/Asymptotics of the Tstat of variance - AC, Homogeneous variance.pdf")

grid.arrange(
  g16,
  g17,
  g18,
  g19,
  g20,
  nrow = 3, ncol = 2,
  top = "Asymptotics of the Tstat of variance - AC, Homogeneous variance", 
  bottom = textGrob(
    "this footnote is right-justified",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  ))
dev.off()
