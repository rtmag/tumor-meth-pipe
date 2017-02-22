
library(RnBeads)
rnb.set.norm=load.rnb.set("rnb.set.norm.RData.zip")

anno=read.table("sample.annotation_tp53_kras_braf.txt",sep="\t",header=T)
rnb.set.norm@pheno=anno

rnb.set.norm@pheno$race=as.character(rnb.set.norm@pheno$tp53_info)
rnb.set.norm@pheno$race[rnb.set.norm@pheno$race!="Normal"]="WT"

norm_dmr <- rnb.execute.computeDiffMeth(rnb.set.norm,pheno.cols=c("race"))

comparison <- get.comparisons(norm_dmr)[1]
norm_dmr_table <-get.table(norm_dmr, comparison, "sites", return.data.frame=TRUE)

meth.norm<-meth(rnb.set.norm)

x=which(norm_dmr_table$diffmeth.p.adj.fdr<0.05 & abs(norm_dmr_table$mean.diff)>.15 )

meth.norm.sig=meth.norm[x,]

tp53=as.character(anno$tp53_info)
tp53[tp53=='Normal']=1
tp53[tp53=='Mutant']=2
tp53[tp53=='WT']=3
tp53=as.numeric(tp53)

braf=as.character(anno$braf_info)
braf[braf=='Normal']=1
braf[braf=='Mutant']=2
braf[braf=='WT']=3
braf=as.numeric(braf)

kras=as.character(anno$kras_info)
kras[kras=='Normal']=1
kras[kras=='Mutant']=2
kras[kras=='WT']=3
kras=as.numeric(kras)

colores=c("white","red","black")
clab=cbind(colores[tp53],colores[braf],colores[kras])
colnames(clab)=c("TP53","BRAF","KRAS")

colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
cols=brewer.pal(3, "Set1")


# set the custom distance and clustering functions
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

#png("COAD_FDR-1_methDif-05.png")
x=heatmap.3(meth.norm.sig,col=colors, hclustfun=hclustfunc, distfun=distfunc, 
            scale="none", trace="none",cexCol=0.2,KeyValueName="Methylation Level",
             ColSideColors=clab,Colv=F,dendrogram="both")
#dev.off()
