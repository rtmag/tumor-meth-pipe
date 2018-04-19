suppressMessages(library(RnBeads))
suppressMessages(library(RColorBrewer))
suppressMessages(library(factoextra))
rnb.set.norm=load.rnb.set("rnb.set.norm.RData.zip")
rnb.set.norm_noNormal=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which(rnb.set.norm@pheno$tp53_info=="Normal")])
noNormal_dmr <- rnb.execute.computeDiffMeth(rnb.set.norm_noNormal,pheno.cols=c("tp53_info"))

comparison <- get.comparisons(noNormal_dmr)[1]
noNormal_dmr_table <-get.table(noNormal_dmr, comparison, "sites", return.data.frame=TRUE)

################# thr FDR diffmeth
meth.norm<-meth(rnb.set.norm)
colnames(meth.norm) = as.character(rnb.set.norm@pheno$tp53_info)
rownames(meth.norm) = rownames(rnb.set.norm@sites)
#meth.mval = mval(rnb.set.norm)
#rownames(meth.mval) = rownames(rnb.set.norm@sites)
#colnames(meth.mval) = as.character(rnb.set.norm@pheno$tp53_info)
meth.norm.sig=meth.norm[noNormal_dmr_table$diffmeth.p.adj.fdr<0.05 & abs(noNormal_dmr_table[,3])>.10,]
saveRDS(meth.norm.sig,"meth.norm.sig.rds")


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
