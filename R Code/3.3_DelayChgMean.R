library(MASS)
library(threejs)
library(amap)
library(igraph)
library(gSeg)

rm(list=ls())

# Import critical value table
CV_table = read.table("/Users/yangwen/Documents/Research/Programing/Data/3.1_OnlineHomoVar_B1000_n.csv",
                      header = TRUE,
                      sep = ",")
B = 1000 # number of test
N = 100 # N is the length of the detection period
M = length(CV_table$d) # M is the numbers of window
nDim = 6 # number of dimensions
nWin = M/nDim # number of windows

# Set up location of change point for detection
change_point = round(N/2)+1 # location of change point
n1 = change_point-1; n2 = N - change_point +1

I_true = rep(1,B)


TPm = 0; FPm = 0; FNm = 0; TNm = 0
TPa = 0; FPa = 0; FNa = 0; TNa = 0

AccuracyMST = rep(0,nDim)
SensitivityMST = rep(0,nDim)
FPR_MST = rep(0,nDim)
Alarm_MST = rep(0,nDim)

AccuracyAC = rep(0,nDim)
SensitivityAC = rep(0,nDim)
FPR_AC = rep(0,nDim)
Alarm_AC = rep(0,nDim)

CV_M1 = rep(0,M)
CV_M2 = rep(0,M)
CV_A1 = rep(0,M)
CV_A2 = rep(0,M)

d = rep(0,nDim)
n = rep(0,nWin)

for (r in 1:nDim) {
  d[r] = CV_table$d[1+(r-1)*nWin]
}
for (r in 1:nWin) {
  n[r] = CV_table$n[r]
} 

for (m in 1:M) {
  CV_M1[m] = CV_table$CU[m]
  CV_M2[m] = CV_table$CU2[m]
  CV_A1[m] = CV_table$ACU[m]
  CV_A2[m] = CV_table$ACU2[m]
}


# Detection procedure 
for (k in 3:nDim) {
  TPm = 0; FPm = 0; FNm = 0; TNm = 0
  TPa = 0; FPa = 0; FNa = 0; TNa = 0
  I_M= rep(0,B)
  I_A= rep(0,B)
  TimeLocation_M= matrix(0,2,B)
  TimeLocation_A= matrix(0,2,B)
  for (i in 1:B){
    # Set multivariate distribution
    mu1 = rep(0, d[k])
    Sigma1 = diag(1,d[k])
    mu2 = rep(0+I_true[i]/((d[k])^(1/3)), d[k])
    Sigma2 = diag(1,d[k]) 
    
    # Generate multivariate distribution
    bivn1 = mvrnorm(n1, mu = mu1, Sigma = Sigma1 ) 
    bivn2 = mvrnorm(n2, mu = mu2, Sigma = Sigma2 )  # from Mass package
    bivn=rbind(bivn1,bivn2)
    
    Record_MST = rep(0,N-n[1]+1)
    Record_AC = rep(0,N-n[1]+1)
    
    # Detection procedure within each window frame: n
    for (j in n[1]:N)
    {
      Win_testM = rep(0,nWin)
      Win_testA = rep(0,nWin)   
      
      for (q in 1:nWin){
        m = (k-1)*nWin + q
        
        if (j < n[q]){ break}else{
          
          # Capture data within the detection frame 
          if (d[k]>1) {
            Bbivn = bivn[(j-n[q]+1):j,]
            Bbivn1 = bivn[(j-n[q]+1):(j-n[q]/2),]
            Bbivn2 = bivn[(j-n[q]/2+1):j,] 
          } else { 
            Bbivn = bivn[(j-n[q]+1):j]
            Bbivn1 = bivn[(j-n[q]+1):(j-n[q]/2)]
            Bbivn2 = bivn[(j-n[q]/2+1):j]   
          }
          
          # Calculate Norm: library(amap)
          W = Dist(Bbivn)
          W1 = Dist(Bbivn1)
          W2 = Dist(Bbivn2)
          
          # Minimum Spanning Tree: library(igraph)
          g = erdos.renyi.game((n[q]), 1)
          g1 = erdos.renyi.game(n[q]/2, 1)
          g2 = erdos.renyi.game(n[q]/2, 1)
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
          
          if ( (TstatM1 > CV_M1[m]) || (TstatM2 > CV_M2[m])) { 
            Win_testM[q] = j-n[q]/2+1
          }
          if ( (TstatA1 > CV_A1[m]) || (TstatA2 > CV_A2[m])) {
            Win_testA[q] = j-n[q]/2+1
          }
        }
      }  
      Record_MST[j] = max(Win_testM)
      Record_AC[j] = max(Win_testA)            
    }
    if (sum(Record_MST) > 0) {I_M[i]=1 #1 for CP # 0 for no CP
    tao = min(which(Record_MST != 0)) #First detected CP
    TimeLocation_M[1,i] = tao  #Record the Time detect CP
    TimeLocation_M[2,i] = Record_MST[tao] #Record the estimated CP
    }
    if (sum(Record_AC) > 0) {I_A[i]=1
    tao = min(which(Record_AC != 0)) #First detected CP
    TimeLocation_A[1,i] = tao  #Record the Time detect CP
    TimeLocation_A[2,i] = Record_AC[tao] #Record the estimated CP
    }
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
  Alarm_MST[k] = sum(TimeLocation_M[1,])/sum(I_M)
  AccuracyMST[k] = (TPm+TNm)/(TPm+FPm+FNm+TNm)
  SensitivityMST[k] = TPm/(TPm+FNm)
  FPR_MST[k]= FPm/(TNm+FPm)
  
  Alarm_AC[k] = sum(TimeLocation_A[1,])/sum(I_A)
  AccuracyAC[k] = (TPa+TNa)/(TPa+FPa+FNa+TNa)
  SensitivityAC[k] = TPa/(TPa+FNa)
  FPR_AC[k]= FPa/(TNa+FPa)
}
Performance = data.frame(d,AccuracyAC,AccuracyMST,SensitivityAC,SensitivityMST,FPR_AC,FPR_MST,Alarm_AC,Alarm_MST)
#saveRDS(ChgPtAC,file="/Users/yangwen/Documents/Research/Programing/Data/ChangePointLocationB100N200AC_17052019.Rda")
write.table(Performance, file = "/Users/yangwen/Documents/Research/Programing/Data/3.3_DelayChgMeanB1000N100_Part2_23052019.csv", sep = ",", col.names = NA,
            qmethod = "double")
