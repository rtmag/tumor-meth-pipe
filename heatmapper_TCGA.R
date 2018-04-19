
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)

files = list.files(pattern="*.meth.norm.sig.rds")

for(i in 16:length(files) ){

meth.norm.sig = readRDS(files[i])
cancer = gsub(".meth.norm.sig.rds","",files[i])

track=colnames(meth.norm.sig)
track[track=="Mutant"]=1
track[track=="WT"]=2
track[track=="Normal"]=3
track=as.numeric(track)

colores=c("#db4e68","#497bd1","#d1c349")
clab=as.character(colores[track])
  
if (dim(meth.norm.sig)[1]<20) next
if (dim(meth.norm.sig)[1]>15000) {
   differences = rowMeans(meth.norm.sig[,track==2])-rowMeans(meth.norm.sig[,track==1])
    meth.norm.sig = tail(meth.norm.sig[order(differences),],15000)
  }

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
tiffname = paste(cancer,".tiff",sep="")
tiff(file = tiffname, width = 3200, height = 3200, units = "px", res = 800)
x = heatmap.2(meth.norm.sig,col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,labCol = "",xlab="", ylab=paste(dim(meth.norm.sig)[1],"CpGs"),key.title="Methylation lvl",ColSideColors=clab,
             keysize=1, key.par = list(cex=0.4),main = cancer)
dev.off()

#pdf("legend_heatmap.pdf")
#plot.new();legend("center",legend=c("TP53Mutant","WT","Normal"),fill=c("#db4e68","#497bd1","#d1c349"),border=NA,bty = "n",cex = 3)
#dev.off()

p53_pattern=names(which((rowMeans(meth.norm.sig[,track==2])-rowMeans(meth.norm.sig[,track==1]))>0))
rdsname = paste(cancer,"_p53pattern.rds",sep="")
saveRDS(p53_pattern,rdsname)
  }


COAD = readRDS("COAD_p53pattern.rds")
LUSC = readRDS("LUSC_p53pattern.rds")
LUAD = readRDS("LUAD_p53pattern.rds")
HNSC = readRDS("HNSC_p53pattern.rds")
MESO = readRDS("MESO_p53pattern.rds")
SARC = readRDS("SARC_p53pattern.rds")
STAD = readRDS("STAD_p53pattern.rds")
PAAD = readRDS("PAAD_p53pattern.rds")
LAML = readRDS("LAML_p53pattern.rds")
UCEC = readRDS("UCEC_p53pattern.rds")
LIHC = readRDS("LIHC_p53pattern.rds")

cgs = c(COAD,LUSC,LUAD,HNSC,MESO,SARC,STAD,PAAD,LAML,UCEC,LIHC)
cgs = unique(cgs)

cpg = matrix(0L, nrow = length(cgs), ncol = 11)
rownames(cpg) = cgs
colnames(cpg) = c("COAD","LUSC","LUAD","HNSC","MESO","SARC","STAD","PAAD","LAML","UCEC","LIHC")
cpg[rownames(cpg) %in% COAD, "COAD"] = 1
cpg[rownames(cpg) %in% LUSC, "LUSC"] = 1
cpg[rownames(cpg) %in% LUAD, "LUAD"] = 1
cpg[rownames(cpg) %in% HNSC, "HNSC"] = 1
cpg[rownames(cpg) %in% MESO, "MESO"] = 1
cpg[rownames(cpg) %in% SARC, "SARC"] = 1
cpg[rownames(cpg) %in% STAD, "STAD"] = 1
cpg[rownames(cpg) %in% PAAD, "PAAD"] = 1
cpg[rownames(cpg) %in% LAML, "LAML"] = 1
cpg[rownames(cpg) %in% UCEC, "UCEC"] = 1
cpg[rownames(cpg) %in% LIHC, "LIHC"] = 1

cpg = cpg[rowSums(cpg)>1,]
x = heatmap.2(cpg,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="CpG")

cuty=cutree(as.hclust(x$rowDendrogram),k=4)
rlab=c("#ffb3ba","#ffdfba","#baffc9","#bae1ff")[cuty]

tiff(file = 'CpG_heatmap.tiff', width = 3200, height = 3200, units = "px", res = 800)
x = heatmap.2(cpg,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="CpG",RowSideColors=rlab)
dev.off()

library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation.table = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

c1 = names(cuty[cuty==1])
c2 = names(cuty[cuty==2])
c3 = names(cuty[cuty==3])
c4 = names(cuty[cuty==4])

anno_c1= annotation.table[rownames(annotation.table) %in% c1,c(1,2)]
#anno_c1[,2] = anno_c1[,2]-1

anno_c2= annotation.table[rownames(annotation.table) %in% c2,c(1,2)]
#anno_c2[,2] = anno_c2[,2]-1

anno_c3= annotation.table[rownames(annotation.table) %in% c3,c(1,2)]
#anno_c3[,2] = anno_c3[,2]-1

anno_c4= annotation.table[rownames(annotation.table) %in% c4,c(1,2)]
#anno_c4[,2] = anno_c4[,2]-1

write.table(anno_c1, "anno_cpg_c1.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(anno_c2, "anno_cpg_c2.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(anno_c3, "anno_cpg_c3.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(anno_c4, "anno_cpg_c4.bed",sep="\t",quote=F,row.names=F,col.names=F)
