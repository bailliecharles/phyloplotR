
library("ape")

if(!require(geiger)){
  install.packages("geiger") #same for 'geiger'
  library(geiger)
}

list.files()

##### pcg 123
pcg123 <- read.nexus(file="mb_123.nex.con.tre" )
pcg123
pcg123$tip.label

new_tips = sort(new_tips, decreasing = FALSE) # trickier if tip labs not ordered
new_tips 

pcg123$tip.label=new_tips

pcg123<-root(pcg123,"ROOT", 
         node = NULL, resolve.root = FALSE, interactive = FALSE)

pcg123 = ladderize(pcg123)

# nodelabels
class(pcg123$node.label)
pp<- as.numeric(pcg123$node.label)
pp<- round(pp,digits=2)

p <- pp ##create a vector, p, that is identical to pp.

p[pp >= 0.99] <- "#1b9e77"
p[pp < 0.99 & pp >= 0.75] <- "#7570b3"
p[pp < 0.75] <- "#d95f02"
p[is.na(pp)] <- "white"

my_pch = c(4, rep(21,pcg123$Nnode)) # X for root node, circles for as many nodes in pcg123
my_pch

plot.phylo(pcg123, use.edge.length=TRUE, direction=
             "r", show.node.label=FALSE, show.tip.label=TRUE, 
           label.offset=0.1, cex=.85, edge.width = 1.5)

#nodelabels(pp, adj = c(1.5, 1.5), frame="none", cex=.5)
#nodelabels(pch=21, cex=1, bg=p[-1])
nodelabels(pch=my_pch, cex=1.5, bg=p)

add.scale.bar(0.1, 2, cex=.8)





