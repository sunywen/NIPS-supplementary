library(MASS)
library(threejs)
library(amap)
#library(easyGgplot2)
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
Sigma2 = 1*diag(1,d) #matrix(1*c(1, 0, 0, 1), 2)
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

# Plot the histogram: library(easyGgplot22)
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

# # Contour plot overlayed on heat map image of results
   bivn.kde = kde2d(bivn[,1], bivn[,2], n = 50)
   head(bivn)
   
   png(file="/Users/yangwen/Desktop/201707 IRTG renewal/Images/2D_contourM5S4.png",width=1200,height=1050,
       bg="transparent",pointsize=36)
   image(bivn.kde)       # from base graphics package
   contour(bivn.kde, add = TRUE)     # from base graphics package
   dev.off()
 # 3D Javascript plot: library(threejs)
 # Unpack data from kde grid format
  x = bivn.kde$x
  y = bivn.kde$y
  z = bivn.kde$z
  # Construct x,y,z coordinates
  xx = rep(x,times=length(y))
  yy = rep(y,each=length(x))
  zz = z; dim(zz) = NULL
 # Set up color range
 ra = ceiling(16 * zz/max(zz))
 col = rainbow(16, 2/3)
  # 3D interactive scatter plot
 # png(file="/Users/yangwen/Desktop/201707 IRTG renewal/Images/ScatterplotM5S4.png",width=1200,height=1050,
 #     bg="transparent",pointsize=36)
 # scatterplot3js(x=xx,y=yy,z=zz,size=0.45,color = col[ra],bg="lightblue")
 # dev.off()

# dd=300
#Time series plot
# x1 = bivn[,1]
# y1 = bivn[,2]
# z1 = 1:length(x1)
# scatterplot3js(x=x1[1:dd],y=y1[1:dd],z=z1[1:dd],size=0.45,color = col[ra],bg="white",xlim=c(-3.5,11),ylim=c(-3.5,11),zlim=c(0,300))
# 


 # Plot MST
 E(g)$width = 1#*E(g)$weight
 E(g)$color = 'grey'
 V(g)$width = 3
 V(g)$color = "orange"
 V(g)$size = 6
 V(g)$label.cex=1.2
 # V(mst)[1:n1]$color = 'orange'
 # V(mst)[n1+1:n]$color = 'pink'
 E(mst)$width = 3
 E(mst)$color = 'blue'
 V(mst)$width = 3
 V(mst)$color = "orange"
 V(mst)$size = 6
 V(mst)$label.cex=1.2

 layout_g  = layout.kamada.kawai(g)
 layout_mst  = layout.kamada.kawai(mst)
 png(file="/Users/yangwen/Desktop/ERG_Transparent.png",width=1200,height=1050,bg='transparent')
 plot(g, layout=layout_g)
 dev.off()
# par(new=F)
 png(file="/Users/yangwen/Desktop/MST_Transparent.png",width=1200,height=1050,bg='transparent')
 plot(mst, layout=layout_mst)
 dev.off()

# par(new=F)
