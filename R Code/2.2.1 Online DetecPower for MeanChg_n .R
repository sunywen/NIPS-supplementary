library(MASS)
library(threejs)
library(amap)
library(igraph)
library(gSeg)

rm(list=ls())

# Import critical value table
CV_table = read.table("/Users/yangwen/Documents/Research/Programing/Data/OnlinePermuCali with HomoVar_B1000_n.csv",
                      header = TRUE,
                      sep = ",")
B = 100 # number of test
N = 100 # N is the length of the detection period
M = length(CV_table$d) # M is the numbers of window

# Set up location of change point for detection
change_point = round(N/2) # location of change point
n1 = change_point; n2 = N - change_point

AccuracyMST = rep(0,M)
SensitivityMST = rep(0,M)
FPR_MST = rep(0,M)

AccuracyAC = rep(0,M)
SensitivityAC = rep(0,M)
FPR_AC = rep(0,M)

for (m in 1:M){
  # Read dimension, window size, and critical value
  d = CV_table$d[m]
  n = CV_table$n[m] 
  
  CV_M1 = CV_table$CU[m]
  CV_M2 = CV_table$CU2[m]
  CV_A1 = CV_table$ACU[m]
  CV_A2 = CV_table$ACU2[m]
  
  I_true = sample(c(0,1), replace=TRUE, size=B)
  I_testM= matrix(0,B,(N-n+1))
  I_testA= matrix(0,B,(N-n+1))
  I_M= rep(0,B)
  I_A= rep(0,B)
  
  TPm = 0; FPm = 0; FNm = 0; TNm = 0
  TPa = 0; FPa = 0; FNa = 0; TNa = 0

  
  # Detection procedure 
  for (i in 1:B){
    # Set multivariate distribution
    mu1 = rep(0, d)
    Sigma1 = diag(1,d) #matrix(c(1, 0, 0, 1), 2)
    mu2 = rep(0+I_true[i]/(d^(1/3)), d)
    Sigma2 = diag(1,d) #diag(1+I_true[i]*1,d) #* (I_true[i]*3+1)
    
    # Generate multivariate distribution
    bivn1 = mvrnorm(n1, mu = mu1, Sigma = Sigma1 ) 
    bivn2 = mvrnorm(n2, mu = mu2, Sigma = Sigma2 )  # from Mass package
    bivn=rbind(bivn1,bivn2)
    

    # Detection procedure within each window frame: n
    for (j in 1:(N-n+1))
    {
      # Capture data within the detection frame 
      if (d>1) {
        Bbivn = bivn[j:(j+n-1),]
        Bbivn1 = bivn[j:(j+n/2-1),]
        Bbivn2 = bivn[(j+n/2):(j+n-1),] 
      } else { 
        Bbivn = bivn[j:(j+n-1)]
        Bbivn1 = bivn[j:(j+n/2-1)]
        Bbivn2 = bivn[(j+n/2):(j+n-1)]   
      }
      
      # Calculate Norm: library(amap)
      W = Dist(Bbivn)
      W1 = Dist(Bbivn1)
      W2 = Dist(Bbivn2)

      # Minimum Spanning Tree: library(igraph)
      g = erdos.renyi.game((n), 1)
      g1 = erdos.renyi.game(n/2, 1)
      g2 = erdos.renyi.game(n/2, 1)
      E(g)$weight = W
      E(g1)$weight = W1
      E(g2)$weight = W2
      mst = minimum.spanning.tree(g)
      mst1 = minimum.spanning.tree(g1)
      mst2 = minimum.spanning.tree(g2)
      
      # Calculate Test statistics
      # MST
      SD = sum(E(mst)$weight)
      SD1 = sum(E(mst1)$weight)
      SD2 = sum(E(mst2)$weight)
      TstatM1 = SD / (SD1 + SD2)
      TstatM2 = SD1 / SD2 + SD2 / SD1
   
      # Complete Graph (AC)
      ASD = sum(E(g)$weight)
      ASD1 = sum(E(g1)$weight)
      ASD2 = sum(E(g2)$weight)
      TstatA1 = ASD / (ASD1 + ASD2)
      TstatA2 = ASD1 / ASD2 + ASD2 / ASD1
  
      if ( (TstatM1 > CV_M1) || (TstatM2 > CV_M2)) { 
          I_testM[i,j] = 1
          }
      if ( (TstatA1 > CV_A1) || (TstatA2 > CV_A2)) {
          I_testA[i,j] = 1
          }
    
    }
    if (sum(I_testM[i,]) > 0) {I_M[i]=1}    
    if (sum(I_testA[i,]) > 0) {I_A[i]=1}  
#    tapply(df$Sales, df$ID, function(a)head(which(a>0),1))
    
    # Calculate testing power
    if (I_true[i]==1 & I_M[i]==1){TPm = TPm + 1}
    if (I_true[i]==1 & I_M[i]==0){FNm = FNm + 1}
    if (I_true[i]==0 & I_M[i]==1){FPm = FPm + 1}
    if (I_true[i]==0 & I_M[i]==0){TNm = TNm + 1}

    if (I_true[i]==1 & I_A[i]==1){TPa = TPa + 1}
    if (I_true[i]==1 & I_A[i]==0){FNa = FNa + 1}
    if (I_true[i]==0 & I_A[i]==1){FPa = FPa + 1}
    if (I_true[i]==0 & I_A[i]==0){TNa = TNa + 1}
  }
  
  AccuracyMST[m] = (TPm+TNm)/(TPm+FPm+FNm+TNm)
  SensitivityMST[m] = TPm/(TPm+FNm)
  FPR_MST[m]= FPm/(TNm+FPm)
  
  AccuracyAC[m] = (TPa+TNa)/(TPa+FPa+FNa+TNa)
  SensitivityAC[m] = TPa/(TPa+FNa)
  FPR_AC[m]= FPa/(TNa+FPa)
 }

Performance = data.frame(d=CV_table$d,n=CV_table$n,AccuracyAC,AccuracyMST,SensitivityAC,SensitivityMST,FPR_AC,FPR_MST)
#saveRDS(Performance,file="/Users/yangwen/Documents/Research/Programing/Data/Performance_2.2_ChgMeanB1000_3rt_10052019.Rda")
write.table(Performance, file = "/Users/yangwen/Documents/Research/Programing/Data/Performance_2.2.1_ChgMeanB100_3rt_n_FPR_20052019.csv", sep = ",", col.names = NA,
            qmethod = "double")
