library(MASS)
library(amap)
library(igraph)
library(gSeg)

rm(list=ls())

# Import critical value table
CV_table = read.table("/Users/yangwen/Documents/Research/Programing/Data/Bootstrap_CV_Homo1_n_a025.csv",
                      header = TRUE,
                      sep = ",")
N = 100 # number of test
M = length(CV_table$d)

AccuracyMST = rep(0,M)
SensitivityMST = rep(0,M)
FPR_MST = rep(0,M)

AccuracyAC = rep(0,M)
SensitivityAC = rep(0,M)
FPR_AC = rep(0,M)

AccuracyR = rep(0,M)
SensitivityR = rep(0,M)
FPR_R = rep(0,M)

for (m in 1:M){
  # Read dimension, window size, and critical value
  d = CV_table$d[m]
  n = CV_table$n[m]
  CV_M1 = CV_table$CU[m]
  CV_M2 = CV_table$CU2[m]
  CV_A1 = CV_table$ACU[m]
  CV_A2 = CV_table$ACU2[m]
  
  change_point = n/2 # location of change point
  n1 = change_point; n2 = n - change_point
  
  TstatM1 = rep(0,N) 
  TstatM2 = rep(0,N)
  
  TstatA1 = rep(0,N)
  TstatA2 = rep(0,N)
  
  r1 = rep(0,N)
  I_true = sample(c(0,1), replace=TRUE, size=N)
  I_testM= rep(0,N)
  I_testA= rep(0,N)
  
  TPm = 0; FPm = 0; FNm = 0; TNm = 0
  TPa = 0; FPa = 0; FNa = 0; TNa = 0
  TPr = 0; FPr = 0; FNr = 0; TNr = 0
  
  # Generate multivariate data: library(MASS)
  for (i in 1:N){
    # Set multivariate distribution
    mu1 = rep(0, d)
    Sigma1 = diag(1,d) #matrix(c(1, 0, 0, 1), 2)
    mu2 = rep(0, d)
    Sigma2 = diag(1+I_true[i],d) #* (I_true[i]*3+1)
    bivn1 = mvrnorm(n1, mu = mu1, Sigma = Sigma1 ) 
    bivn2 =  mvrnorm(n2, mu = mu2, Sigma = Sigma2 )  # from Mass package
    bivn=rbind(bivn1,bivn2)
    
    # Calculate Norm: library(amap)
    W = Dist(bivn)
    W1 = Dist(bivn1)
    W2 = Dist(bivn2)
    #norm(as.matrix(bivn[1,]-bivn[2,]),"F")
    
    # Minimum Spanning Tree: library(igraph)
    g = erdos.renyi.game(n, 1)
    g1 = erdos.renyi.game(n1, 1)
    g2 = erdos.renyi.game(n2, 1)
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
    TstatM1[i] = SD / (SD1 + SD2)
    TstatM2[i] =  (SD1 / SD2 + SD2 / SD1)
    
    # All connected (AC)
    ASD = sum(E(g)$weight)
    ASD1 = sum(E(g1)$weight)
    ASD2 = sum(E(g2)$weight)
    TstatA1[i] = ASD / (ASD1 + ASD2)
    TstatA2[i] = ASD1 / ASD2 + ASD2 / ASD1
    
    if ( (TstatM1[i] > CV_M1) || (TstatM2[i] > CV_M2)) { I_testM[i] = 1}
    if ( (TstatA1[i] > CV_A1) || (TstatA2[i] > CV_A2)) { I_testA[i] = 1}
    
    # Calculat CPD by Chen & Zhang
    E = as_edgelist(mst)
    if (gseg1(n, E)$pval.appr < 0.05) { r1[i] = 1 }
    
    # Calculate testing power
    if (I_true[i]==1 & I_testM[i]==1){TPm = TPm + 1}
    if (I_true[i]==1 & I_testM[i]==0){FNm = FNm + 1}
    if (I_true[i]==0 & I_testM[i]==1){FPm = FPm + 1}
    if (I_true[i]==0 & I_testM[i]==0){TNm = TNm + 1}
    
    if (I_true[i]==1 & I_testA[i]==1){TPa = TPa + 1}
    if (I_true[i]==1 & I_testA[i]==0){FNa = FNa + 1}
    if (I_true[i]==0 & I_testA[i]==1){FPa = FPa + 1}
    if (I_true[i]==0 & I_testA[i]==0){TNa = TNa + 1}
    
    if (I_true[i]==1 & r1[i]==1){TPr = TPr + 1}
    if (I_true[i]==1 & r1[i]==0){FNr = FNr + 1}
    if (I_true[i]==0 & r1[i]==1){FPr = FPr + 1}
    if (I_true[i]==0 & r1[i]==0){TNr = TNr + 1}
  }
  
  AccuracyMST[m] = (TPm+TNm)/(TPm+FPm+FNm+TNm)
  SensitivityMST[m] = TPm/(TPm+FNm)
  FPR_MST[m]= FPm/(TNm+FPm)
  
  AccuracyAC[m] = (TPa+TNa)/(TPa+FPa+FNa+TNa)
  SensitivityAC[m] = TPa/(TPa+FNa)
  FPR_AC[m]= FPa/(TNa+FPa)
  
  AccuracyR[m] = (TPr+TNr)/(TPr+FPr+FNr+TNr)
  SensitivityR[m] = TPr/(TPr+FNr)
  FPR_R[m]= FPr/(TNr+FPr)
}

Performance = data.frame(d=CV_table$d,n=CV_table$n,AccuracyAC,AccuracyMST,AccuracyR,SensitivityAC,SensitivityMST,SensitivityR,FPR_AC,FPR_MST,FPR_R)
#saveRDS(Performance,file="/Users/yangwen/Documents/Research/Programing/Data/Performance_Var1p5_10052019.Rda")
write.table(Performance, file = "/Users/yangwen/Documents/Research/Programing/Data/Performance_Var2_n_23052019.csv", sep = ",", col.names = NA,
            qmethod = "double")

cat(AccuracyMST)
cat(SensitivityMST)
cat(AccuracyAC)
cat(SensitivityAC)
cat(AccuracyR)
cat(SensitivityR)

