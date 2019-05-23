library(rasterImage)
library(raster)
library(png)
library(grid)

for (i in 1:3){
  img <- readPNG(paste("/Users/yangwen/Desktop/Contour",i,".png",sep=""))
  png(file=paste("/Users/yangwen/Desktop/distr-",i,".png",sep=""),width=448,height=366)
  grid.raster(img[1:448,1:366,])
  dev.off()
}