suppressMessages(library(RnBeads))
suppressMessages(library(RColorBrewer))

rnb.set.norm=load.rnb.set("rnb.set.norm.RData.zip")
rnb.set.norm_noNormal=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which(rnb.set.norm@pheno$tp53_info=="Normal")])
noNormal_dmr <- rnb.execute.computeDiffMeth(rnb.set.norm_noNormal,pheno.cols=c("tp53_info"))

comparison <- get.comparisons(noNormal_dmr)[1]
noNormal_dmr_table <-get.table(noNormal_dmr, comparison, "sites", return.data.frame=TRUE)

################# thr FDR diffmeth
meth.norm<-meth(rnb.set.norm)
meth.norm.sig=meth.norm[noNormal_dmr_table$diffmeth.p.adj.fdr<0.01 & abs(noNormal_dmr_table[,3])>.05,]
colnames(meth.norm.sig)<-rnb.set.norm@pheno$tp53_info
colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
cols=brewer.pal(3, "Set1")

# set the custom distance and clustering functions
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

# perform clustering on rows and columns (genes and samples)
cl.row <- hclustfunc(distfunc(meth.norm.sig))
gr.row <- cutree(cl.row, 2)

cl.col <- hclustfunc(distfunc(t(meth.norm.sig)))
colOrder <- cl.col$labels
colOrder[colOrder=="Mutant"]=1
colOrder[colOrder=="WT"]=2
colOrder[colOrder=="Normal"]=3
colOrder=as.numeric(colOrder)

col1 <- c("brown","orange")

tiff(file = "TP53_wtVSmt_FDR-1_methDif-05.tiff", width = 3200, height = 3200, units = "px", res = 800) 
x=heatmap.2(meth.norm.sig,col=colors, hclustfun=hclustfunc, distfun=distfunc, scale="none", trace="none",cexCol=0.2,ColSideColors=cols[colOrder],RowSideColors=col1[gr.ro
w],lhei=c(2,4), lwid=c(2,3.5), keysize=0.75, key.par = list(cex=0.3))
dev.off()
