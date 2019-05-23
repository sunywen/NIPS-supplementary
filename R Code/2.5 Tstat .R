library(MASS)
library(threejs)
library(amap)
library(igraph)
library(gSeg)
library(ggplot2)

rm(list=ls())

# Import critical value table
CV_table = read.table("/Users/yangwen/Documents/Research/Programing/Data/OnlinePermuCali with HomoVar_B1000_n.csv",
                      header = TRUE,
                      sep = ",")
B = 100 # number of test
N = 200 # N is the length of the detection period
M = length(CV_table$d) # M is the numbers of window

# Set up location of change point for detection
change_point = round(N/2)+1 # location of change point
n1 = change_point-1; n2 = N - change_point +1

SensitivityMST = rep(0,M)
ChgPtMST = rep(0,B)
ChgPtErMST = rep(0,B)
DelayMST = rep(0,M)
SensitivityAC = rep(0,M)
ChgPtAC = rep(0,B)
ChgPtErAC = rep(0,B)
DelayAC = rep(0,M)


for (m in 40:40){
  # Read dimension, window size, and critical value
  d = CV_table$d[m]
  n = CV_table$n[m] 
  
  CV_M1 = CV_table$CU[m]
  CV_M2 = CV_table$CU2[m]
  CV_A1 = CV_table$ACU[m]
  CV_A2 = CV_table$ACU2[m]
  
  I_M= rep(0,B)
  I_A= rep(0,B)
  
  # Detection procedure 
  for (i in 1:B){
    
    I_testM = rep(0,(N-n+1))
    I_testA = rep(0,(N-n+1))
    TstatM1 = rep(0,(N-n+1))
    TstatM2 = rep(0,(N-n+1))
    TstatA1 = rep(0,(N-n+1))
    TstatA2 = rep(0,(N-n+1))
    
    # Set multivariate distribution
    mu1 = rep(0, d)
    Sigma1 = diag(1,d) #matrix(c(1, 0, 0, 1), 2)
    mu2 = rep(0, d) #rep(1/(d^(1/3)), d)
    Sigma2 = diag(2,d) #diag(1+I_true[i]*1,d) #* (I_true[i]*3+1)
    
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
      TstatM1[j] = SD / (SD1 + SD2)
      TstatM2[j] = SD1 / SD2 + SD2 / SD1
      
      # Complete Graph (AC)
      ASD = sum(E(g)$weight)
      ASD1 = sum(E(g1)$weight)
      ASD2 = sum(E(g2)$weight)
      TstatA1[j] = ASD / (ASD1 + ASD2)
      TstatA2[j] = ASD1 / ASD2 + ASD2 / ASD1
      
      if ( (TstatM1[j] > CV_M1) || (TstatM2[j] > CV_M2)) { 
        I_testM[j] = (j+n/2)
      }
      if ( (TstatA1[j] > CV_A1) || (TstatA2[j] > CV_A2)) {
        I_testA[j] = (j+n/2)
      }
      
    }
    if (sum(I_testM) > 0) {I_M[i]=1
                           ChgPtMST[i] = min(I_testM[I_testM > 0])
                           ChgPtErMST[i] = abs(ChgPtMST[i]-change_point)
    } # Output change point location  
    if (sum(I_testA) > 0) {I_A[i]=1
                           ChgPtAC[i] = min(I_testA[I_testA > 0])
                           ChgPtErAC[i] = abs(ChgPtAC[i]-change_point)
    }  # Output change point location
                                                                
    
  
  }
  
  SensitivityMST[m] = sum(I_M)/B
  SensitivityAC[m] = sum(I_A)/B
  DelayMST[m] = sum(ChgPtErMST)/sum(I_M)
  DelayAC[m] = sum(ChgPtErAC)/sum(I_A)
  
}
Tstamp = (1+n/2):(N-n/2+1)
TestStat = data.frame(Tstamp, TstatA1,TstatA2,TstatM1,TstatM2)

#Plot

ggplot(data = TestStat, aes(x = Tstamp, y = TstatA2))+
  geom_line(color = "#00AFBB", size = 2)

cat("ChgPtMST", ChgPtMST,"\n")
cat("ChgPtAC", ChgPtAC,"\n")
cat("Sensitivity MST:",SensitivityMST[m]," Sensitivity AC:",SensitivityAC[m],"\n")
cat("DelayMST:", DelayMST[m], "Delay AC:", DelayAC[m],"\n")

#saveRDS(TestStat,file="/Users/yangwen/Documents/Research/Programing/Data/Performance_2.2_ChgMeanB1000_3rt_10052019.Rda")
#write.table(Performance, file = "/Users/yangwen/Documents/Research/Programing/Data/2.5 Online test statistic along time_d100_n30.csv", sep = ",", col.names = NA,
#            qmethod = "double")
