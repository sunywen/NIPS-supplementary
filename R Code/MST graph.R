library(igraph)
g = erdos.renyi.game(20, 0.8)
E(g)$weight=runif(length(E(g)), 0., 1.)
E(g)$width= 3.0*E(g)$weight

mst = minimum.spanning.tree(g)
E(mst)$color = 'blue'

layout  = layout.kamada.kawai(mst)

plot(g, layout=layout)
par(new=T)
E(mst)$width  =  3
plot(mst, layout=layout)
par(new=F)

