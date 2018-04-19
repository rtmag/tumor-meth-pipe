args = commandArgs(trailingOnly=TRUE)
suppressMessages(library(RnBeads))
filename=paste(args[1],"rnb.set.norm.RData.zip",sep="")
rnb.set.norm=load.rnb.set(filename)
rnb.set.norm_noNormal=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which(rnb.set.norm@pheno$tp53_info=="Normal")])
noNormal_dmr <- rnb.execute.computeDiffMeth(rnb.set.norm_noNormal,pheno.cols=c("tp53_info"))

comparison <- get.comparisons(noNormal_dmr)[1]
noNormal_dmr_table <-get.table(noNormal_dmr, comparison, "sites", return.data.frame=TRUE)

################# thr FDR diffmeth
meth.norm<-meth(rnb.set.norm)
colnames(meth.norm) = as.character(rnb.set.norm@pheno$tp53_info)
rownames(meth.norm) = rownames(rnb.set.norm@sites)
meth.norm.sig=meth.norm[noNormal_dmr_table$diffmeth.p.adj.fdr<0.05 & abs(noNormal_dmr_table[,3])>.10,]

name=gsub("/root/TCGA/TCGA/preprocessing_tp53/","",args[1])
name=gsub("/","",name)
name=paste("/root/TCGA/TCGA/matrix/",name,".meth.norm.sig.rds",sep="")

saveRDS(meth.norm.sig,name)
#



