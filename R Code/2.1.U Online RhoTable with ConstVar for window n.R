library(MASS)
library(threejs)
library(amap)
library(igraph)
library(gSeg)

rm(list=ls())

# Import critical value table
CV_table = read.table("/Users/yangwen/Documents/Research/Programing/Data/Critical_value_n.csv",
                      header = TRUE,
                      sep = ",")
B = 1000 # B is the number of bootstrap
M = length(CV_table$d)
alpha = 0.05

CL=rep(0,M)
CU=rep(0,M)


ACL=rep(0,M)
ACU=rep(0,M)


for (m in 1:M){
  # Read dimension, window size, and critical value
  d = CV_table$d[m]
  n = CV_table$n[m]  
  
  N = 100 # N is the length of the training sample
  
  change_point = round(N/2) # location of change point
  n1 = change_point; n2 = N - change_point
  
  # Set multivariate distribution
  mu1 = rep(0, d)
  Sigma1 = diag(1,d) 
  mu2 = rep(0, d)
  Sigma2 = diag(1,d) 
  Tstat = rep(0,B)
  Tdiff = rep(0,B)

  TstatA = rep(0,B)
  TdiffA = rep(0,B)
 
  # Generate multivariate data: library(MASS)
  bivn1 = mvrnorm(n1, mu = mu1, Sigma = Sigma1) 
  bivn2 =  mvrnorm(n2, mu = mu2, Sigma = Sigma2)  
  bivn=rbind(bivn1,bivn2)
  
  
  # Permutation procedure 
  for (i in 1:B){
    RdS = sample(1:nrow(bivn), replace = FALSE)
    # Note: replace= FALSE performs permutation and it is bootstrap in our case, which finds the MST from different combination of sample nodes
    
    if (d>1) {
      BbivnRds = bivn[RdS,]
    } else {
      BbivnRds = bivn[RdS]
    }
    
    tTstat = rep(0,(N-n+1))

    tTstatA = rep(0,(N-n+1))
 
    for (j in 1:(N-n+1))
    {
      
      if (d>1) {
        Bbivn = BbivnRds[j:(j+n-1),]
        Bbivn1 = BbivnRds[j:(j+n/2-1),]
        Bbivn2 = BbivnRds[(j+n/2):(j+n-1),] 
      } else { 
        Bbivn = BbivnRds[j:(j+n-1)]
        Bbivn1 = BbivnRds[j:(j+n/2-1)]
        Bbivn2 = BbivnRds[(j+n/2):(j+n-1)]   
      }
      
      # Calculate Norm: library(amap)
      
      g = erdos.renyi.game((n), 1)
      W = Dist(Bbivn)
      E(g)$weight = W
      mst = minimum.spanning.tree(g)
      SD = sum(E(mst)$weight)
      ASD = sum(E(g)$weight)
      
      W1 = Dist(Bbivn1)
      W2 = Dist(Bbivn2)
      
      # Minimum Spanning Tree: library(igraph)
      g1 = erdos.renyi.game(n/2, 1)
      g2 = erdos.renyi.game(n/2, 1)
      E(g1)$weight = W1
      E(g2)$weight = W2
      mst1 = minimum.spanning.tree(g1)
      mst2 = minimum.spanning.tree(g2)
      
      # Calculate Test statistics
      # MST
      SD1 = sum(E(mst1)$weight)
      SD2 = sum(E(mst2)$weight)
      tTstat[j] = SD / (SD1 + SD2)*(SD1 / SD2 + SD2 / SD1)
      
      
      # All connected
      ASD1 = sum(E(g1)$weight)
      ASD2 = sum(E(g2)$weight)
      tTstatA[j] = ASD / (ASD1 + ASD2)*(ASD1 / ASD2 + ASD2 / ASD1)

    }
    Tstat[i] = max(tTstat)

    TstatA[i] = max(tTstatA)

    
  }
  
  # Calculate the quantile at (alpha, 1-alpha)
  # MST
  Tstat_hat = mean(Tstat)
  Tdiff=Tstat-Tstat_hat
  seTb= sqrt(1/(B-1)*sum(Tdiff^2))
  CL[m] = Tstat_hat + seTb * quantile(Tdiff/seTb, alpha)
  CU[m] = Tstat_hat + seTb * quantile(Tdiff/seTb, 1-alpha)
  cat(" T1 Quantiles at ", alpha,1-alpha, "are ", CL[m], CU[m])
  

  # All connected
  Tstat_hatA = mean(TstatA)
  TdiffA=TstatA-Tstat_hatA
  seTbA= sqrt(1/(B-1)*sum(TdiffA^2))
  ACL[m] = Tstat_hatA + seTbA * quantile(TdiffA/seTbA, alpha)
  ACU[m] = Tstat_hatA + seTbA * quantile(TdiffA/seTbA, 1-alpha)
  cat(" T1A Quantiles at ", alpha,1-alpha, "are ", ACL[m], ACU[m])
  
  

  # Plot the histogram of bootstrap statistics
  #ggplot2.histogram(data=Tstat, xName='weight',
  #                 fill="white", color="black",
  #                  addDensityCurve=TRUE, densityFill='#FF6666')
}

Performance = data.frame(d=CV_table$d,n=CV_table$n,CL,CU,ACL,ACU)
#saveRDS(Performance,file="/Users/yangwen/Desktop/Research/CV_table.Rda")
write.table(Performance, file = "/Users/yangwen/Documents/Research/Programing/Data/OnlinePermuCali with HomoVar_B1000_n_UTstat.csv", sep = ",", col.names = NA,
            qmethod = "double")
