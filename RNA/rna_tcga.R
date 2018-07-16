library(gplots)
library(ggplot2)
library(graphics)
library(scales)
options(scipen=999)
library(gplots)
library(RColorBrewer)
library("xlsx")

rna_samples = read.table("sample.txt",sep='\t')
rna_samples = cbind(rna_samples,rna_samples)

sample_450k = read.table("450k_annotation_tcga_p53.txt",sep='\t',header=TRUE)

#match(a,b) # a[ !is.na(match) ]  # b[which(match)]
ix= match(rna_samples[,1], sample_450k[,1])

#head(rna_samples[!is.na(ix),])
#head(sample_450k[ix[!is.na(ix)],])

rna = rownames(rna_samples)[!is.na(ix)]
k450 = sample_450k[ix[!is.na(ix)],]



rpkm <- do.call(cbind, lapply(rna, function(x) read.table(pipe(paste('zcat',x) ), sep="\t", row.names=1 ) ))

ix_foxm1 = grep("ENSG00000111206",rownames(rpkm))

FOXM1 = data.frame(cell = as.character(k450[,2]), gene=log2(as.numeric(rpkm[ix_foxm1,])) )
pdf("TCGA_FOXM1.pdf")
boxplot(gene ~ cell, vertical = TRUE, data = FOXM1,ylab = expression('Log2 FOXM1 FPKM'),names=c("TP53 Mutant","Normal","TP53 WT"),outline=F)
stripchart(gene ~ cell, vertical = TRUE, data = FOXM1, jitter = 0.3, ylab = expression('Log2 FOXM1 FPKM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2,add=T)
dev.off()
#ENSG00000124762
ix_CDKN1 = grep("ENSG00000124762",rownames(rpkm))
                              
CDKN1 = data.frame(cell = as.character(k450[,2]), gene=log2(as.numeric(rpkm[ix_CDKN1,])) )
pdf("TCGA_CDKN1.pdf")
boxplot(gene ~ cell, vertical = TRUE, data = CDKN1,ylab = expression('Log2 CDKN1A FPKM'),names=c("TP53 Mutant","Normal","TP53 WT"),outline=F)
stripchart(gene ~ cell, vertical = TRUE, data = CDKN1, jitter = 0.3, ylab = expression('Log2 CDKN1A FPKM'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.5),cex = 2,add=T)
dev.off()

####################################################################################
mut = read.csv("pheno_pChange.csv")
id = t(matrix(unlist(unlist(strsplit(as.character(mut$Sample_ID),"-"))),nrow=7))
id = id[,1:4]
id =  apply(id, 1, function(x) paste(x, collapse="-"))
id = cbind(id, as.character(mut$Protein_Change))
            
##########################
mtNO =  as.numeric(gsub("\\D", "", id[,2]))
id[which(mtNO>99 & mtNO < 290),2] = "DNA Binding Domain"
id[which(mtNO>318),2] = "Tetra motif"
id[grep("p.",id[,2]),2] = "Other Mutations"
id = id[match(k450[,1],id[,1]),]
##########################
id[grep("p.R175",id[,2]),2] = "R175"
id[grep("p.R273",id[,2]),2] = "R273"
id[grep("p.R248",id[,2]),2] = "R248"
id[grep("p.",id[,2]),2] = "Other Mutations"
id = id[match(k450[,1],id[,1]),]
            
FOXM1 = data.frame( cell = as.character(k450[,2]), gene=log2(as.numeric(rpkm[ix_foxm1,])), pChange =id[,2]  )

                              
                              
boxplot(gene ~ cell, vertical = TRUE, data = FOXM1,ylab = expression('Log2 CDKN1A FPKM'),names=c("TP53 Mutant","Normal","TP53 WT"),outline=F)
stripchart(gene ~ cell, vertical = TRUE, data = FOXM1, jitter = 0.3, ylab = expression('Log2 CDKN1A FPKM'),
    method = "jitter", pch = 20, col = c(rep('red',200),rep('blue',23)),cex = 2,add=T)

pdf("FOXM1_domains.pdf")
ggplot(FOXM1, aes(x=pChange, y=gene, color=pChange)) +
  geom_jitter(position=position_jitter(0.2))
dev.off()

pdf("FOXM1_hotspot.pdf")
ggplot(FOXM1, aes(x=pChange, y=gene, color=pChange)) +
  geom_jitter(position=position_jitter(0.2))
dev.off()
                              
