library(RnBeads)
rnb.set.norm=load.rnb.set("rnb.set.norm.RData.zip")

anno=read.table("sample.annotation_tp53_kras_braf.txt",sep="\t",header=T)
rnb.set.norm@pheno=anno

rnb.set.norm_noNormal=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which(rnb.set.norm@pheno$tp53_info=="Normal")])

tp53_dmr <- rnb.execute.computeDiffMeth(rnb.set.norm_noNormal,pheno.cols=c("tp53_info"))
braf_dmr <- rnb.execute.computeDiffMeth(rnb.set.norm_noNormal,pheno.cols=c("braf_info"))
kras_dmr <- rnb.execute.computeDiffMeth(rnb.set.norm_noNormal,pheno.cols=c("kras_info"))

comparison <- get.comparisons(tp53_dmr)[1]
tp53_dmr_table <-get.table(tp53_dmr, comparison, "sites", return.data.frame=TRUE)

comparison <- get.comparisons(braf_dmr)[1]
braf_dmr_table <-get.table(braf_dmr, comparison, "sites", return.data.frame=TRUE)

comparison <- get.comparisons(kras_dmr)[1]
kras_dmr_table <-get.table(kras_dmr, comparison, "sites", return.data.frame=TRUE)

meth.norm<-meth(rnb.set.norm)

x=which(  tp53_dmr_table$diffmeth.p.adj.fdr<0.01 & abs(tp53_dmr_table$mean.diff)>.05 )

meth.norm.sig=meth.norm[x,]


tp53=as.character(anno$tp53_info)
tp53[tp53=='Normal']=1
tp53[tp53=='Mutant']=2
tp53[tp53=='WT']=3
tp53=as.numeric(tp53)



colores=c("white","red","black")
clab=cbind(colores[tp53])
colnames(clab)=c("TP53")

colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
cols=brewer.pal(3, "Set1")


# set the custom distance and clustering functions
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

#tiff("COAD_FDR-1_methDif-05.tiff",res = 300)
colnames(meth.norm.sig)=NULL
rownames(meth.norm.sig)=NULL
x=heatmap.3(meth.norm.sig,col=colors, hclustfun=hclustfunc, distfun=distfunc, 
            scale="none", trace="none",cexCol=0.2,KeyValueName="Methylation Level",
             ColSideColors=clab,Colv=T,dendrogram="none",lhei=c(1,15))
            
legend("topright",legend=c("Normal","Mutant","Wild Type"),
fill=c(colores), border=T, bty="n", y.intersp = 0.7, cex=0.7)
           
dev.off()
