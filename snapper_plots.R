source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")

library(ggtree)
library(ape)
library(geiger)


getwd()
setwd("~/Dropbox/scripts")
list.files()

# Messing about with command line to remove []&label= from figtree file for import with node labels
# tr -d Œ[]&=,' < in.file | sed 's/label//g' > out.file

# Read in tree
twee = read.newick("twee_nwk_MASTER") #rooted tree from figtree. 
plot(twee, show.tip.label = T)
twee

# ladderize with root at bottom
twee = ladderize(twee, right=FALSE)

# proportional branch lengths
twee2 = compute.brlen(twee) # twee2 becomes master
# get rid of annoying tips
twee2 = drop.tip(twee2, "Tip1") 
twee2 = drop.tip(twee2,"Tip2")

plot(twee2, show.tip.label = F)

twee2$tip.label

#### APE METHOD ####

#### tip labs
cols  = read.csv("colours.csv", header=F) # file has: V1=families, V2=taxon names, V3=hex codes for colours of groups
class(cols)
head(cols)
cols = cols[,1:3]
head(cols)
cols$V2 <- gsub(">", "", cols$V2)  #remove arrow
cols = cols[order(cols$V2),]       #reorder alphbetical
row.names(cols)=c(1:nrow(cols))    #reorder

# sort rows of df to match order of tip labs
# this seems really convoluted but spec chars causing issues with ordering????? Need to come up with better solution. 
twee2$tip.label <- gsub("[[:punct:]]", "", twee2$tip.label) 
x = order(twee2$tip.label) # order
cols3 = reorder(cols$V2, x)
cols4 = order(cols3)
cols2 = cols[order(x),] # reorder df based on index
cols2
write.csv(cols2, "colours_orderedNames.csv")# write out new df 
write.tree(twee2, file="TEMP_twee2_newNames") # write out temp tree to check tiplabs

# Give temp tip labels
twee2$tip.label = cols2$V2


### colors for clades

# nodes 

X1 = tips(twee2, node=338) # 338 
X2 = tips(twee2, node=341)#341 
X3 = tips(twee2, node=331)#331
X4= tips(twee2, node=332) #332
X5 = tips(twee2, node=423)#423
X6= tips(twee2, node=297) #297
X7 = tips(twee2, node=327)#327
X8= tips(twee2, node=326)#326
X9 = tips(twee2, node=323)#323
X10 = tips(twee2, node=316)#316


X1branches = which.edge(twee2, X1)
X2branches = which.edge(twee2, X2)
X3branches = which.edge(twee2, X3)
X4branches= which.edge(twee2, X4)
X5branches = c(which.edge(twee2, X5),347) # adding other branch colouring
X6branches= which.edge(twee2, X6)
X7branches = which.edge(twee2, X7)
X8branches = which.edge(twee2, X8)
X9branches =which.edge(twee2, X9)
X10branches = which.edge(twee2, X10)

clcolr <- rep("#4F5F69", dim(twee2$edge)[1])
clcolr[X1] = "#245790"
clcolr[X2] = "#793d8c"
clcolr[X3] = "#C43541"
clcolr[X4] = "#8C234D"
clcolr[X5] = "#2B2373"
clcolr[X6] = "#608b26"
clcolr[X7] = "#8f6b6b"
clcolr[X8] = "#621e05"
clcolr[X9] = "#ab4326"
clcolr[X10] = "#f39940"

# The above can be done more neatly by putting everything in df and making a for loop


plot(twee2, edge.color=clcolr, edge.width=2, show.tip.label=F)

### node labs
a<- as.numeric(twee2$node.label)
            
a[a<=69] <-0
a[a>=70] <- "X"
a
nodelabels(pch=a, cex=0.45, col="black")
twee2$node.label = a

###= names
cols3 = read.csv("colours_orderedNames.csv", header=T) # had to move out to parse italic part: FUCKING awkward names
	# has: V1=full names except those with species where just accession and '|', V2=species names, blank for else, V3= '|' and location for those with sp names, else blank. 

#new tip labels 
twee2$tip.label <- mixedFontLabel(as.character(cols3$V1), as.character(cols3$V2), as.character(cols3$V3), italic = 2)
twee2$tip.label

#### FINAL Plotting

plot(twee2, cex=0.45, label.offset=0.01, edge.color=clcolr, edge.width=1.3, show.tip.label=T, no.margin=T)
nodelabels(pch=a, cex=0.4, bg=p, col="black")
add.scale.bar(length=0.1, lwd=0.4, cex=0.4, ask=T )
title(main="what do you want for a title?", cex.main = 1.5, col.main="pink", line = -7, adj=0.1)
title(main="and/or a subtitle?", cex.main = 1, col.main="red", line = -8, adj=0.1)



############################
##### GGTREE Version #######

# first sort out names 
cols3 = read.csv("colours_orderedNames.csv", header=T) # name file as above
cols3 = cols3[,2:4]
labels = twee2$tip.label
cols3 = data.frame(label=labels,cols3)
apply(cols3, 2, class)
cols3$V2 <- paste0("'", cols3$V2, "'") # need to surround contents in quotes for parsing tip labs later
cols3$V3 <- paste0("'", cols3$V3, "'")
cols3$V4 <- paste0("'", cols3$V4, "'")


twee3 = groupClade(twee2, node=c(338,341,331,332,423,297,327,326,323,316,311,312,256,255,254,250)) # group clades

myCols = c("#4F5F69","#245790", "#793d8c" ,"#C43541" ,"#8C234D", "#2B2373" ,	# colors for branches
           "#608b26", "#8f6b6b", "#621e05", "#ab4326" ,"#f39940", "#c0d818",
           "#c99615" ,"#068054", "#38c19a","#799689", "#98a5ad")

twee4 = ggtree(twee3, color="#4F5F69", size=1.5) %<+% cols3 + xlim(0,1.4) + 	# plot, attach data, and size						
		geom_tiplab(aes(label=paste0('plain(', V2, ')~italic(', V3, ')~', V4)), parse=T, cex=2.5, offset=0.03, color="black")+ #names
			scale_color_manual(values=myCols) # color branches
# The above plots a tree with grey branches. To coor the branches by clade change color="#45F69" to  aes(color=group)

# can also add bars by family 
twee4+
geom_cladelabel(node=338, label=" ",  align=T,barsize=2, color="#245790", fontsize = 2.5)+
geom_cladelabel(node=341, label=" ",  align=T, barsize=2, color="#793d8c", fontsize = 2.5)+
geom_cladelabel(node=331, label=" ",  align=T, barsize=2, color="#C43541", fontsize = 2.5)+
geom_cladelabel(node=332, label=" ",  align=T, barsize=2, color="#8C234D", fontsize = 2.5)+
geom_cladelabel(node=423, label=" ",  align=T, barsize=2, color="#2B2373", fontsize = 2.5)+
geom_cladelabel(node=297, label=" ", align=T, barsize=2, color="#608b26", fontsize = 2.5)+
geom_cladelabel(node=327, label=" ", align=T, barsize=2, color="#8f6b6b", fontsize =2.5)+
geom_cladelabel(node=326, label=" ", align=T, barsize=2, color="#621e05", fontsize = 2.5)+
geom_cladelabel(node=323, label=" ", align=T,  barsize=2, color="#ab4326", fontsize = 2.5)+
geom_cladelabel(node=316, label=" ", align=T,  barsize=2, color="#f39940", fontsize = 2.5)+
# add a scale
geom_treescale(x=0.4,y=30, color="#4F5F69", fontsize = 3, offset=1.5)+
# add node labels - only for those above 0.7 support
	# need to make this a symbol... 
geom_text2(aes(label=bootstrap, subset=bootstrap>0.7))

bootstrap=twee2$node.label











                     
                     
                     

