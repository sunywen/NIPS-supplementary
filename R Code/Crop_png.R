library(rasterImage)
library(raster)
library(png)
library(grid)

for (i in 0:32){
img <- readPNG(paste("/Users/yangwen/Desktop/Distributions/",i,".png",sep=""))
png(file=paste("/Users/yangwen/Desktop/Distributions/distr-",i,".png",sep=""),width=900,height=625)
grid.raster(img[50:1300,400:2200,])
dev.off()
}