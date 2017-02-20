 rnb.set.norm=load.rnb.set("rnb.set.norm.zip")
meth.norm<-meth(rnb.set.norm)
colnames(meth.norm)=(rnb.set.norm@pheno$Sample_ID)
meth.norm2=meth.norm[,1:16]
sig=readRDS("~/tblab/TCGA/preprocessing_tp53/COAD/FDR1p_5pmeth_diff.rds")
gr.row=readRDS("~/tblab/TCGA/preprocessing_tp53/COAD/gr.row.rds")

meth.norm.sig=meth.norm2[sig,]
colnames(meth.norm.sig)<-colnames(meth.norm2)
colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
cols=brewer.pal(3, "Set1")

# set the custom distance and clustering functions
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

# perform clustering on rows and columns (genes and samples)
cl.row <- hclustfunc(distfunc(meth.norm.sig))
gr.row <- cutree(cl.row, 2)

col1 <- c("brown","orange")

#png("COAD_FDR-1_methDif-05.png")
x=heatmap.2(meth.norm.sig,col=colors, hclustfun=hclustfunc, distfun=distfunc, scale="none", trace="none",cexCol=0.2,RowSideColors=col1[gr.row],Colv=F,dendrogram="row")
#dev.off()

png("COAD_DiffCpG_450Stephanie.png")
x=heatmap.2(meth.norm.sig,col=colors, hclustfun=hclustfunc, distfun=distfunc, scale="none", trace="none",cexCol=.7,Colv=F,dendrogram="row")
dev.off()


png("COAD_hypo_DiffCpG_450Stephanie.png")
x=heatmap.2(meth.norm.sig[gr.row==2,],col=colors, hclustfun=hclustfunc, distfun=distfunc, scale="none", trace="none",cexCol=.7,Colv=F,dendrogram="row")
dev.off()
