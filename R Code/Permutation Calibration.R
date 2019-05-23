# Goal: find critical value for online CPD
# Definition: training sample size: N, fixed window size: n, dimension: d, 
#             and permutation steps: B
# Note: Calibration is performed for all window frame captured within the trainning sample
#       yet only sigle window size n is applied    

library(MASS)
library(threejs)
library(amap)
library(easyGgplot2)
library(igraph)
library(gSeg)

rm(list=ls())

# Import critical value table
CV_table = read.table("/Users/yangwen/Desktop/Research/Critical_value.csv",
                      header = TRUE,
                      sep = ",")

B = 10000 # B is the number of permutation
M = length(CV_table$d)
alpha = 0.05

CL=rep(0,M)
CU=rep(0,M)

for (m in 1:M){
  # Read dimension, window size, and critical value
  d = CV_table$d[m]
  n = CV_table$n[m]  
  
  # Set multivariate distribution
  mu1 = rep(0, d)
  Sigma1 = diag(1,d) 
  
  Tstat = rep(0,B)
  Tdiff = rep(0,B)
  
  # Generate multivariate data: library(MASS)
  bivn = mvrnorm(n1, mu = mu1, Sigma = Sigma1) 

  # window size n
  # for 
  W = Dist(bivn)
  g = erdos.renyi.game(n, 1)
  E(g)$weight = W
  mst = minimum.spanning.tree(g)
  SD = sum(E(mst)$weight)
  
  # Permutation procedure 
  for (i in 1:B){
    RdS = sample(1:nrow(bivn), replace = FALSE)
    #Note: replace= FALSE performs permutation, which finds the MST from different combination of sample nodes
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
    
    # Minimum Spanning Tree: library(igraph)
    g1 = erdos.renyi.game(n1, 1)
    g2 = erdos.renyi.game(n2, 1)
    E(g1)$weight = W1
    E(g2)$weight = W2
    mst1 = minimum.spanning.tree(g1)
    mst2 = minimum.spanning.tree(g2)
    
    # Calculate Test statistics
    SD1 = sum(E(mst1)$weight)
    SD2 = sum(E(mst2)$weight)
    Tstat[i] = SD * (SD1 / SD2 + SD2 / SD1) / (SD1 + SD2)
  }
  Tstat_hat = mean(Tstat)
  Tdiff=Tstat-Tstat_hat
  seTb= sqrt(1/(B-1)*sum(Tdiff^2))
  CL[m] = Tstat_hat + seTb * quantile(Tdiff/seTb, alpha)
  CU[m] = Tstat_hat + seTb * quantile(Tdiff/seTb, 1-alpha)
  cat("Quantiles at ", alpha,1-alpha, "are ", CL[m], CU[m])
  
  # Plot the histogram of bootstrap statistics  
  #ggplot2.histogram(data=Tstat, xName='weight',
  #                 fill="white", color="black",
  #                  addDensityCurve=TRUE, densityFill='#FF6666')
}

CV_tableN = data.frame(d=CV_table$d,n=CV_table$n,CL,CU)
saveRDS(CV_tableN,file="/Users/yangwen/Desktop/Research/CV_table.Rda")
write.table(Performance, file = "/Users/yangwen/Desktop/Research/Bootstrap_CV.csv", sep = ",", col.names = NA,
            qmethod = "double")
