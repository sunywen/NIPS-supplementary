library(MASS)
library(threejs)
library(amap)
library(easyGgplot2)
library(igraph)
library(gSeg)

rm(list=ls())

# Set change point
N = 1 #number of simulations
n = 60; #window size
d = 2
change_point = n/2 # location of change point
n1 = change_point; n2 = n - change_point

# Set multivariate distribution
mu1 = rep(0, d)
Sigma1 = diag(1,d) #matrix(c(1, 0, 0, 1), 2)
mu2 = rep(0, d)
Sigma2 = 4*diag(1,d) #matrix(1*c(1, 0, 0, 1), 2)
#cat("Multivariate data: mu1 = ", mu1, "\n")
#cat("                   sigma1 = ", Sigma1, "\n")
#cat("Multivariate data: mu2 = ", mu2, "\n")
#cat("                   sigma2 = ", Sigma2, "\n")
Tstat = rep(0,N)
TstatA= rep(0,N)
r1 = rep(0,N)

# Generate multivariate data: library(MASS)
for (i in 1:N){
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
  SD = sum(E(mst)$weight)
  SD1 = sum(E(mst1)$weight)
  SD2 = sum(E(mst2)$weight)
  Tstat[i] = SD / (SD1 + SD2) 

  # Calculat CPD by Chen & Zhang
  E = as_edgelist(mst)
  r1[i] = gseg1(n, E)$Zmax
}

# Plot the histogram: library(easyGgplot2)
#ggplot2.histogram(data=TstatA, xName='weight',
#                  fill="white", color="black",
#                  addDensityCurve=TRUE, densityFill='#FF6666')

# ggplot2.histogram(data=Tstat, xName='weight',
#                   fill="white", color="black",
#                   addDensityCurve=TRUE, densityFill='#FF6666')

# ggplot2.histogram(data=r1, xName='weight',
#                  fill="white", color="black",
#                  addDensityCurve=TRUE, densityFill='#FF6666')

Tq = quantile(Tstat, c(0.05, 0.95))
Rq = quantile(r1, c(0.05, 0.95))
cat("Critical value for Tstat at 5%:", Tq,"\n") 
cat("Critical value for r1 at 5%:", Rq,"\n")

# # # Contour plot overlayed on heat map image of results
#    bivn.kde = kde2d(bivn[,1], bivn[,2], n = 50)
#    head(bivn)
#    image(bivn.kde)       # from base graphics package
#    contour(bivn.kde, add = TRUE)     # from base graphics package
# # # 3D Javascript plot: library(threejs)
# # # Unpack data from kde grid format
#  x = bivn.kde$x
#  y = bivn.kde$y
#  z = bivn.kde$z
#  # Construct x,y,z coordinates
#  xx = rep(x,times=length(y))
#  yy = rep(y,each=length(x))
#  zz = z; dim(zz) = NULL
# # # Set up color range
#  ra = ceiling(16 * zz/max(zz))
#  col = rainbow(16, 2/3)
# # # 3D interactive scatter plot
#  scatterplot3js(x=xx,y=yy,z=zz,size=0.4,color = col[ra],bg="black")
#  
#  
# #Time series plot
#  x1 = bivn[,1]
#  y1 = bivn[,2]
#  z1 = 1:length(x1)
#  scatterplot3js(x=x1[1:30],y=y1[1:30],z=z1[1:30],size=0.4,color = col[ra],bg="black")
# # 

 # Plot MST
 E(g)$width = 0.5*E(g)$weight
 V(g)$size = 6
 E(mst1)$width=0.5
 E(mst1)$label.cex=1
 E(mst1)$color = 'blue'
 E(mst1)$width = 3
 V(mst1)$color = 'orange'
 V(mst1)$width = 1
 V(mst1)$size = 6
 V(mst1)$label.cex=1
 
 E(mst2)$width=0.5
 E(mst2)$label.cex=1
 E(mst2)$color = 'green'
 E(mst2)$width = 3
 V(mst2)$color = 'red'
 V(mst2)$width = 1
 V(mst2)$size = 6
 V(mst2)$label.cex=1
 
 layout1  = layout.kamada.kawai(mst1)
 layout2  = layout.kamada.kawai(mst2)
 #png(file="/Users/yangwen/Desktop/Research Presentation/Images/G3.png",width=1200,height=1050)
 plot(mst2, layout=layout2)
 #dev.off()
 #par(new=F)
 #png(file="/Users/yangwen/Desktop/Research Presentation/Images/MST_CPsigma3.png",width=1200,height=1050)
 plot(mst1, layout=layout1)
 #dev.off()
 
# par(new=T)
