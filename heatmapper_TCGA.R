
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)
meth.norm.sig = readRDS("meth.norm.sig.rds")

track=colnames(meth.norm.sig)
track[track=="Mutant"]=1
track[track=="WT"]=2
track[track=="Normal"]=3
track=as.numeric(track)

colores=c("#db4e68","#497bd1","#d1c349")
clab=as.character(colores[track])

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
x = heatmap.2(meth.norm.sig,col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,labCol = "",xlab="", ylab="CpGs",key.title="Methylation lvl",ColSideColors=clab)
