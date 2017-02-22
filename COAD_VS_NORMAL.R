
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

#x=which(norm_dmr_table$diffmeth.p.adj.fdr<0.05 & abs(norm_dmr_table$mean.diff)>.25 )

x=order(norm_dmr_table$combinedRank)[1:30000]


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
##

# stage
stage=as.character(anno$ajcc_pathologic_tumor_stage)
stage[grep("Stage IIA",stage)]="Stage II"
stage[grep("Stage IIB",stage)]="Stage II"
stage[grep("Stage IIC",stage)]="Stage II"
stage[grep("Stage IVB",stage)]="Stage IV"
stage[grep("Stage IIIA",stage)]="Stage III"
col.stage=brewer.pal(8, "Set1")
col.stage[1]="white"
stage.names=names(table(stage))

stage[stage=='[Not Available]']=1
stage[stage=='Stage I']=2
stage[stage=='Stage II']=3
stage[stage=='Stage III']=4
stage[stage=='Stage IIIB']=5
stage[stage=='Stage IIIC']=6
stage[stage=='Stage IV']=7
stage[stage=='Stage IVA']=8

stage=as.numeric(stage)
##

clab=cbind(colores[tp53],colores[braf],colores[kras],col.stage[stage])
colnames(clab)=c("TP53","BRAF","KRAS","Stage")


###

colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)

# set the custom distance and clustering functions
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

#png("COAD_FDR-1_methDif-05.png")
x=heatmap.3(meth.norm.sig,col=colors, hclustfun=hclustfunc, distfun=distfunc, 
            scale="none", trace="none",cexCol=0.2,KeyValueName="Methylation Level",
             ColSideColors=clab,dendrogram="both")
            
legend("topright",legend=c("Normal","Mutant","Wild Type",stage.names),
fill=c(colores,col.stage), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
#dev.off()
