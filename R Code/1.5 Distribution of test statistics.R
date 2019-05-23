library(MASS)
library(threejs)
library(amap)
library(igraph)
library(ggplot2)
library(gSeg)
library(grid)
library(gridExtra)

rm(list=ls())

# Import critical value table
CV_table = read.table("/Users/yangwen/Desktop/Research/Programing/Data/Critical_value.csv",
                      header = TRUE,
                      sep = ",")

B = 1000 # B is the number of runs
M = length(CV_table$d)
alpha = 0.05

ACL=rep(0,M)
ACU=rep(0,M)
ACL2=rep(0,M)
ACU2=rep(0,M)

TstatA = matrix(0,M,B)
TstatA2 = matrix(0,M,B)

for (m in 1:M){
  # Read dimension, window size, and critical value
  d = CV_table$d[m]
  n = CV_table$n[m]  
  change_point = n/2 # location of change point
  n1 = change_point; n2 = n - change_point
  
  # Set multivariate distribution
  mu1 = rep(0, d)
  Sigma1 = diag(1,d) 
  mu2 = rep(0, d)
  Sigma2 = diag(1,d)
  

  # Generate multivariate data: library(MASS)
  bivn1 = mvrnorm(n1, mu = mu1, Sigma = Sigma1) 
  bivn2 =  mvrnorm(n2, mu = mu2, Sigma = Sigma2)  
  bivn=rbind(bivn1,bivn2)
  
  W = Dist(bivn)
  g = erdos.renyi.game(n, 1)
  E(g)$weight = W
  ASD = sum(E(g)$weight)
  
  # Permutation procedure 
  for (i in 1:B){
    RdS = sample(1:nrow(bivn), replace = FALSE)
    #Note: replace= FALSE performs permutation and it is bootstrap in our case, which finds the MST from different combination of sample nodes
    if (d>1) {
      Bbivn = bivn[RdS,]
      Bbivn1 = Bbivn[1:n1,]
      Bbivn2 = Bbivn[n1+1:n2,] 
    } else { 
      Bbivn = bivn[RdS]
      Bbivn1 = Bbivn[1:n1]
      Bbivn2 = Bbivn[n1+1:n2]   
    }
    
    # Calculate Norm: library(amap)
    W1 = Dist(Bbivn1)
    W2 = Dist(Bbivn2)
    g1 = erdos.renyi.game(n1, 1)
    g2 = erdos.renyi.game(n2, 1)
    E(g1)$weight = W1
    E(g2)$weight = W2

    # Calculate Test statistics
   
    # All connected
    ASD1 = sum(E(g1)$weight)
    ASD2 = sum(E(g2)$weight)
    TstatA[m,i] = ASD / (ASD1 + ASD2)
    TstatA2[m,i] = ASD1 / ASD2 + ASD2 / ASD1
  }
 
  # All connected
  
  ACL[m] = quantile(TstatA[m,], alpha)
  ACU[m] = quantile(TstatA[m,], 1-alpha)
  cat(" T1A Quantiles at ", alpha,1-alpha, "are ", ACL[m], ACU[m])

  ACL2[m] = quantile(TstatA2, alpha)
  ACU2[m] = quantile(TstatA2, 1-alpha)
  cat(" T2A Quantiles at ", alpha,1-alpha, "are ", ACL2[m], ACU2[m])

}

#Plot the histogram of test statistics
# Histogram for the "AGE" column in the "chol" dataset, with title "Histogram for Age" and label for the x-axis ("Age"), with bins of a width of 0.5 that range from values 20 to 50 on the x-axis and that have transparent blue filling and red borders
qplot(TstatA[3,],
      geom="histogram",
      #binwidth = 0.1,  
      main = "Histogram for Age", 
      xlab = "Age",  
      fill=I("blue"), 
      col=I("red"), 
      alpha=I(.2),
     # xlim=c(2.0,2.1)
      )

# Plot emperical distributions
Tm1 = data.frame(stat=TstatA[1,])
Tm2 = data.frame(stat=TstatA[2,])
Tm3 = data.frame(stat=TstatA[3,])
Tm4 = data.frame(stat=TstatA[4,])
Tm5 = data.frame(stat=TstatA[5,])
Tm6 = data.frame(stat=TstatA[6,])
Tm7 = data.frame(stat=TstatA[7,])
Tm8 = data.frame(stat=TstatA[8,])
Tm9 = data.frame(stat=TstatA[9,])
Tm10 = data.frame(stat=TstatA[10,])
Tm11 = data.frame(stat=TstatA[11,])
Tm12 = data.frame(stat=TstatA[12,])
Tm13 = data.frame(stat=TstatA[13,])
Tm14 = data.frame(stat=TstatA[14,]) 
Tm15 = data.frame(stat=TstatA[15,])

Tm1$d = 'd=1'
Tm2$d = 'd=1'
Tm3$d = 'd=1'
Tm4$d = 'd=10'
Tm5$d = 'd=10'
Tm6$d = 'd=10'
Tm7$d = 'd=50'
Tm8$d = 'd=50'
Tm9$d = 'd=50'
Tm10$d = 'd=100'
Tm11$d = 'd=100'
Tm12$d = 'd=100'
Tm13$d = 'd=500'
Tm14$d = 'd=500'
Tm15$d = 'd=500'

Tm1$n = 'n=15'
Tm2$n = 'n=30'
Tm3$n = 'n=100'
Tm4$n = 'n=15'
Tm5$n = 'n=30'
Tm6$n = 'n=100'
Tm7$n = 'n=15'
Tm8$n = 'n=30'
Tm9$n = 'n=100'
Tm10$n = 'n=15'
Tm11$n = 'n=30'
Tm12$n = 'n=100'
Tm13$n = 'n=15'
Tm14$n = 'n=30'
Tm15$n = 'n=100'

Tv1 = data.frame(stat=TstatA2[1,])
Tv2 = data.frame(stat=TstatA2[2,])
Tv3 = data.frame(stat=TstatA2[3,])
Tv4 = data.frame(stat=TstatA2[4,])
Tv5 = data.frame(stat=TstatA2[5,])
Tv6 = data.frame(stat=TstatA2[6,])
Tv7 = data.frame(stat=TstatA2[7,])
Tv8 = data.frame(stat=TstatA2[8,])
Tv9 = data.frame(stat=TstatA2[9,])
Tv10 = data.frame(stat=TstatA2[10,])
Tv11 = data.frame(stat=TstatA2[11,])
Tv12 = data.frame(stat=TstatA2[12,])
Tv13 = data.frame(stat=TstatA2[13,])
Tv14 = data.frame(stat=TstatA2[14,]) 
Tv15 = data.frame(stat=TstatA2[15,])

Tv1$d = 'd=1'
Tv2$d = 'd=1'
Tv3$d = 'd=1'
Tv4$d = 'd=10'
Tv5$d = 'd=10'
Tv6$d = 'd=10'
Tv7$d = 'd=50'
Tv8$d = 'd=50'
Tv9$d = 'd=50'
Tv10$d = 'd=100'
Tv11$d = 'd=100'
Tv12$d = 'd=100'
Tv13$d = 'd=500'
Tv14$d = 'd=500'
Tv15$d = 'd=500'

Tv1$n = 'n=15'
Tv2$n = 'n=30'
Tv3$n = 'n=100'
Tv4$n = 'n=15'
Tv5$n = 'n=30'
Tv6$n = 'n=100'
Tv7$n = 'n=15'
Tv8$n = 'n=30'
Tv9$n = 'n=100'
Tv10$n = 'n=15'
Tv11$n = 'n=30'
Tv12$n = 'n=100'
Tv13$n = 'n=15'
Tv14$n = 'n=30'
Tv15$n = 'n=100'

# and combine into your new data frame
Mean_d1 =  rbind(Tm1, Tm2, Tm3)
pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of mean, d=1.pdf")
ggplot(Mean_d1, aes(stat, fill = n)) + geom_density(alpha = 0.2) + labs(x="T1", y="density",
       title="Emperical distribution of Test statistics of mean, d=1")
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of mean, d=10.pdf")
Mean_d10 =  rbind(Tm4, Tm5, Tm6)
ggplot(Mean_d10, aes(stat, fill = n)) + geom_density(alpha = 0.2) + labs(x="T1", y="density",
                                                                        title="Emperical distribution of test statistics of mean, d=10")
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of mean, d=50.pdf")
Mean_d50 =  rbind(Tm7, Tm8, Tm9)
ggplot(Mean_d50, aes(stat, fill = n)) + geom_density(alpha = 0.2) + labs(x="T1", y="density",
                                                                        title="Emperical distribution of test statistics of mean, d=50")
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of mean, d=100.pdf")
Mean_d100 =  rbind(Tm10, Tm11, Tm12)
ggplot(Mean_d100, aes(stat, fill = n)) + geom_density(alpha = 0.2) + labs(x="T1", y="density",
                                                                         title="Emperical distribution of test statistics of mean, d=100")
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of mean, d=500.pdf")
Mean_d500 =  rbind(Tm13, Tm14, Tm15)
ggplot(Mean_d500, aes(stat, fill = n)) + geom_density(alpha = 0.2) + labs(x="T1", y="density",
                                                                         title="Emperical distribution of test statistics of mean, d=500")
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of mean, n=15.pdf")
Mean_n15 = rbind(Tm1,Tm4,Tm7,Tm10,Tm13)
ggplot(Mean_n15, aes(stat, fill = d)) + geom_density(alpha = 0.2) + labs(x="T1", y="density",
                                                                          title="Emperical distribution of test statistics of mean, n=15")
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of mean, n=30.pdf")
Mean_n30 = rbind(Tm2,Tm5,Tm8,Tm11,Tm14)
ggplot(Mean_n30, aes(stat, fill = d)) + geom_density(alpha = 0.2) + labs(x="T1", y="density",
                                                                         title="Emperical distribution of test statistics of mean, n=30")
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of mean, n=100.pdf")
Mean_n100 = rbind(Tm3,Tm6,Tm9,Tm12,Tm15)
ggplot(Mean_n100, aes(stat, fill = d)) + geom_density(alpha = 0.2) + labs(x="T1", y="density",
                                                                         title="Emperical distribution of test statistics of mean, n=100")
dev.off()

# Test statistics for variance
pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of variance, d=1.pdf")
Var_d1 =  rbind(Tv1, Tv2, Tv3)
ggplot(Var_d1, aes(stat, fill = n)) + geom_density(alpha = 0.2) + labs(x="T2", y="density",
                                                                        title="Emperical distribution of Test statistics of variance, d=1")
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of variance, d=10.pdf")
Var_d10 =  rbind(Tv4, Tv5, Tv6)
ggplot(Var_d10, aes(stat, fill = n)) + geom_density(alpha = 0.2) + labs(x="T2", y="density",
                                                                         title="Emperical distribution of test statistics of variance, d=10")
dev.off()


pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of variance, d=50.pdf")
Var_d50 =  rbind(Tv7, Tv8, Tv9)
ggplot(Var_d50, aes(stat, fill = n)) + geom_density(alpha = 0.2) + labs(x="T2", y="density",
                                                                         title="Emperical distribution of test statistics of variance, d=50")
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of variance, d=100.pdf")
Var_d100 =  rbind(Tv10, Tv11, Tv12)
ggplot(Var_d100, aes(stat, fill = n)) + geom_density(alpha = 0.2) + labs(x="T2", y="density",
                                                                          title="Emperical distribution of test statistics of variance, d=100")
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of variance, d=500.pdf")
Var_d500 =  rbind(Tv13, Tv14, Tv15)
ggplot(Var_d500, aes(stat, fill = n)) + geom_density(alpha = 0.2) + labs(x="T2", y="density",
                                                                          title="Emperical distribution of test statistics of variance, d=500")
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of variance, n=15.pdf")
Var_n15 = rbind(Tv1,Tv4,Tv7,Tv10,Tv13)
ggplot(Var_n15, aes(stat, fill = d)) + geom_density(alpha = 0.2) + xlim(NA,2.0025) + labs(x="T2", y="density",
                                                                         title="Emperical distribution of test statistics of variance, n=15")
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of variance, n=30.pdf")
Var_n30 = rbind(Tv2,Tv5,Tv8,Tv11,Tv14)
ggplot(Var_n30, aes(stat, fill = d)) + geom_density(alpha = 0.2) + xlim(NA,2.0025) + labs(x="T2", y="density",
                                                                         title="Emperical distribution of test statistics of variance, n=30")
dev.off()

pdf("/Users/yangwen/Desktop/Research/Analysis/Emperical distribution of Test statistics of variance, n=100.pdf")
Var_n100 = rbind(Tv3,Tv6,Tv9,Tv12,Tv15)
ggplot(Var_n100, aes(stat, fill = d)) + geom_density(alpha = 0.2) + xlim(NA,2.0025) + labs(x="T2", y="density",
                                                                          title="Emperical distribution of test statistics of variance, n=100")
dev.off()

# comparison plots
g1=ggplot(Mean_d10, aes(stat, fill = n)) + geom_density(alpha = 0.2) + labs(x="T1", y="density",
                                                     title="Emperical distribution of test statistics of mean, d=10")
g2=ggplot(Var_d10, aes(stat, fill = n)) + geom_density(alpha = 0.2) + labs(x="T2", y="density",
                                                     title="Emperical distribution of test statistics of variance, d=10")

pdf("/Users/yangwen/Desktop/Research/Analysis/Asymptotics of the Tstat of mean and variance.pdf")

grid.arrange(
  g1,
  g2,
  nrow = 2, ncol = 1,
  #top = "Asymptotics of the Tstat - MST and Complete graph, Homogeneous variance",
  bottom = textGrob(
    "",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  ))
dev.off()


g3=ggplot(Mean_n100, aes(stat, fill = d)) + geom_density(alpha = 0.2) + labs(x="T1", y="density",
                                                                          title="Emperical distribution of test statistics of mean, n=100")
g4=ggplot(Var_n100, aes(stat, fill = d)) + geom_density(alpha = 0.2) + xlim(NA,2.0025) + labs(x="T2", y="density",
                                                                                           title="Emperical distribution of test statistics of variance, n=100")
pdf("/Users/yangwen/Desktop/Research/Analysis/Asymptotics of the Tstat of mean and variance- various with d.pdf")

grid.arrange(
  g3,
  g4,
  nrow = 2, ncol = 1,
  #top = "Asymptotics of the Tstat - MST and Complete graph, Homogeneous variance",
  bottom = textGrob(
    "",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 1,
    x = 1
  ))
dev.off()


Performance = data.frame(d=CV_table$d,n=CV_table$n/2,ACL,ACU,ACL2,ACU2)
saveRDS(Performance,file="/Users/yangwen/Desktop/Research/Programing/Data/DistribtuionTstat.Rda")
write.table(Performance, file = "/Users/yangwen/Desktop/Research/Programing/Data/DistributionTstat.csv", sep = ",", col.names = NA,
            qmethod = "double")
