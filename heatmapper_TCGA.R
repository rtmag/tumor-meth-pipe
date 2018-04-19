
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)

files = list.files(pattern="*.meth.norm.sig.rds")

for(i in 1:length(files) ){

meth.norm.sig = readRDS(files[i])
cancer = gsub(".meth.norm.sig.rds","",files[i])

if (dim(meth.norm.sig)[1]<20) next
  
track=colnames(meth.norm.sig)
track[track=="Mutant"]=1
track[track=="WT"]=2
track[track=="Normal"]=3
track=as.numeric(track)

colores=c("#db4e68","#497bd1","#d1c349")
clab=as.character(colores[track])

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
