#Higher Dimensional Distributions
library(corrplot)
library(clusterGeneration)
mu = rep(0,10)
pdMat = genPositiveDefMat(10,lambdaLow=10)
Sigma = pdMat$Sigma
dim(Sigma)
mvn = mvrnorm(5000, mu = mu, Sigma = Sigma )

corrplot(cor(mvn), 
         method="ellipse",
         tl.pos="n",
         title="Matrix Correlations")
