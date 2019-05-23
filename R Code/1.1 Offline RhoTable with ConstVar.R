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

B = 10000 # B is the number of bootstrap
M = length(CV_table$d)
alpha = 0.05

CL=rep(0,M)
CU=rep(0,M)
CL2=rep(0,M)
CU2=rep(0,M)

ACL=rep(0,M)
ACU=rep(0,M)
ACL2=rep(0,M)
ACU2=rep(0,M)

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
  Tstat = rep(0,B)
  Tdiff = rep(0,B)
  Tstat2 = rep(0,B)
  Tdiff2 = rep(0,B)
  
  TstatA = rep(0,B)
  TdiffA = rep(0,B)
  TstatA2 = rep(0,B)
  TdiffA2 = rep(0,B)
  
  # Generate multivariate data: library(MASS)
  bivn1 = mvrnorm(n1, mu = mu1, Sigma = Sigma1) 
  bivn2 =  mvrnorm(n2, mu = mu2, Sigma = Sigma2)  
  bivn=rbind(bivn1,bivn2)
  
  W = Dist(bivn)
  g = erdos.renyi.game(n, 1)
  E(g)$weight = W
  mst = minimum.spanning.tree(g)
  SD = sum(E(mst)$weight)
  ASD = sum(E(g)$weight)
  # Bootstrap procedure 
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
    
    # Minimum Spanning Tree: library(igraph)
    g1 = erdos.renyi.game(n1, 1)
    g2 = erdos.renyi.game(n2, 1)
    E(g1)$weight = W1
    E(g2)$weight = W2
    mst1 = minimum.spanning.tree(g1)
    mst2 = minimum.spanning.tree(g2)
    
    # Calculate Test statistics
    # MST
    SD1 = sum(E(mst1)$weight)
    SD2 = sum(E(mst2)$weight)
    Tstat[i] = SD / (SD1 + SD2)
    Tstat2[i] = SD1 / SD2 + SD2 / SD1
    
    
    # All connected
    ASD1 = sum(E(g1)$weight)
    ASD2 = sum(E(g2)$weight)
    TstatA[i] = ASD / (ASD1 + ASD2)
    TstatA2[i] = ASD1 / ASD2 + ASD2 / ASD1
      }
  # MST
  Tstat_hat = mean(Tstat)
  Tdiff=Tstat-Tstat_hat
  seTb= sqrt(1/(B-1)*sum(Tdiff^2))
  CL[m] = Tstat_hat + seTb * quantile(Tdiff/seTb, alpha)
  CU[m] = Tstat_hat + seTb * quantile(Tdiff/seTb, 1-alpha)
  cat(" T1 Quantiles at ", alpha,1-alpha, "are ", CL[m], CU[m])
  
  Tstat_hat2 = mean(Tstat2)
  Tdiff2=Tstat2-Tstat_hat2
  seTb2= sqrt(1/(B-1)*sum(Tdiff2^2))
  CL2[m] = Tstat_hat2 + seTb2 * quantile(Tdiff2/seTb2, alpha)
  CU2[m] = Tstat_hat2 + seTb2 * quantile(Tdiff2/seTb2, 1-alpha)
  cat(" T2 Quantiles at ", alpha,1-alpha, "are ", CL2[m], CU2[m])
  
  # All connected

  Tstat_hatA = mean(TstatA)
  TdiffA=TstatA-Tstat_hatA
  seTbA= sqrt(1/(B-1)*sum(TdiffA^2))
  ACL[m] = Tstat_hatA + seTbA * quantile(TdiffA/seTbA, alpha)
  ACU[m] = Tstat_hatA + seTbA * quantile(TdiffA/seTbA, 1-alpha)
  cat(" T1A Quantiles at ", alpha,1-alpha, "are ", ACL[m], ACU[m])
  
   
  Tstat_hatA2 = mean(TstatA2)
  TdiffA2=TstatA2-Tstat_hatA2
  seTbA2= sqrt(1/(B-1)*sum(TdiffA2^2))
  ACL2[m] = Tstat_hatA2 + seTbA2 * quantile(TdiffA2/seTbA2, alpha)
  ACU2[m] = Tstat_hatA2 + seTbA2 * quantile(TdiffA2/seTbA2, 1-alpha)
  cat(" T2A Quantiles at ", alpha,1-alpha, "are ", ACL2[m], ACU2[m])
  
  
  # Plot the histogram of bootstrap statistics  
  #ggplot2.histogram(data=Tstat, xName='weight',
  #                 fill="white", color="black",
  #                  addDensityCurve=TRUE, densityFill='#FF6666')
}

Performance = data.frame(d=CV_table$d,n=CV_table$n,CL,CU, CL2,CU2,ACL,ACU,ACL2,ACU2)
saveRDS(Performance,file="/Users/yangwen/Desktop/Research/CV_table.Rda")
write.table(Performance, file = "/Users/yangwen/Desktop/Research/Bootstrap_CV_Homo1.csv", sep = ",", col.names = NA,
            qmethod = "double")
