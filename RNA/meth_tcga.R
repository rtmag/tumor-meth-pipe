options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)
library(RnBeads)

rnb.set.norm=load.rnb.set("rnb.set.norm.RData.zip")
 
probe.annotation <- rnb.get.annotation("probes450")

x=as.data.frame(unlist(unlist(probe.annotation)))

rownames(x) = gsub("chr\\d\\.","",rownames(x),perl =T)

#chr12: 2986099 2986480


colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
pdf("TEST.pdf")
sig_vsd = meth.norm[which(x$seqnames=="chr12" & x$start>2985000 & x$start<2987000),]
colnames(sig_vsd) = (rnb.set.norm@pheno$tp53_info)
  colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
  heatmap.2(sig_vsd,col=colors,scale="none", trace="none",srtCol=90,
labRow = FALSE,xlab="", ylab="Genes",key.title="Methylation",cexCol=.8)
dev.off()
