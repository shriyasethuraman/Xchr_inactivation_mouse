a<- read.table("count_rep1_129",header=F,sep="\t",stringsAsFactors=F)
b<- read.table("count_rep2_129",header=F,sep="\t",stringsAsFactors=F)
c<- read.table("count_rep3_129",header=F,sep="\t",stringsAsFactors=F)
d<- read.table("count_rep4_129",header=F,sep="\t",stringsAsFactors=F)
e<- read.table("count_rep5_single_129",header=F,sep="\t",stringsAsFactors=F)
f<- read.table("count_rep6_single_129",header=F,sep="\t",stringsAsFactors=F)
g<- read.table("count_rep7_single_129",header=F,sep="\t",stringsAsFactors=F)
h<- read.table("count_S12_single_129",header=F,sep="\t",stringsAsFactors=F)
i<- read.table("count_S13_single_129",header=F,sep="\t",stringsAsFactors=F)
j<- read.table("count_S14_single_129",header=F,sep="\t",stringsAsFactors=F)


colnames(a) <- c("gene","WT_129")
colnames(b) <- c("gene","Mut_EZH2")
colnames(c) <- c("gene","Mut_EZH1")
colnames(d) <- c("gene","Mut_EZH1_EZH2")
colnames(e) <- c("gene","Mut_EED_het_129")
colnames(f) <- c("gene","Mut1_EED_129")
colnames(g) <- c("gene","Mut2_EED_129")
colnames(h) <- c("gene","Mut2_EED_het_129")
colnames(i) <- c("gene","Mut3_EED_129")
colnames(j) <- c("gene","Mut4_EED_129")

merge1 <- merge(a,b,by=("gene"))
merge2 <- merge(merge1,c,by=("gene"))
merge3 <- merge(merge2,d,by=("gene"))
merge4 <- merge(merge3,e,by=("gene"))
merge5 <- merge(merge4,f,by=("gene"))
merge6 <- merge(merge5,g,by=("gene"))
merge7 <- merge(merge6,h,by=("gene"))
merge8 <- merge(merge7,i,by=("gene"))
total_count_single<- merge(merge8,j,by=("gene"))

write.table(total_count_single, "total_count_single.txt", sep="\t", row.names=FALSE)

library("DESeq")
CountTable<-read.table ("total_count_single.txt", header=TRUE, row.names=1)
CountTable<-read.table ("total_count_single.txt", header=TRUE, row.names=1)
head(CountTable)

conds <-factor(c("WT","ezh2","ezh1","ezh2_ezh1","eed_het","eed","eed","eed_het","eed","eed"))
conds
cds <- newCountDataSet(CountTable, conds)
head(counts(cds))

#normalization and visualization
cds <- estimateSizeFactors (cds)

sizeFactors(cds)

cds <- estimateDispersions (cds, method="blind")

pdf("genome_dispersions_single.pdf")
(plotDispEsts(cds))
dev.off()

##from here down, repeat for each mutant

res = nbinomTest (cds, "WT", "ezh2")

pdf("plotMA_EZH2_Mut_single.pdf")
print(plotMA(res))
dev.off()

pdf("histogram_EZH2_Mut_single.pdf")
print(hist(res$padj, breaks=100))
dev.off()

resSig <- res[res$padj <0.05,]
head (resSig[order(resSig$pval),])
write.csv(resSig, file="genome_differential_expression_EZH2_Mut_single.csv")
write.csv(res, file="genome_assayed_genes_EZH2_Mut_single.csv")

cdsBlind <- estimateDispersions (cds, method="blind")
vsd <- getVarianceStabilizedData (cdsBlind)
select <-order(res$padj)[1:50]
colors<- colorRampPalette(c("white", "darkblue"))(100)
pdf("genome_heatmap_EZH2_Mut_single.pdf")
print(heatmap(vsd[select,], col=colors, scale="none"))
dev.off()

res = nbinomTest (cds, "WT", "ezh1")

pdf("plotMA_EZH1_Mut_single.pdf")
print(plotMA(res))
dev.off()

pdf("histogram_EZH1_Mut_single.pdf")
print(hist(res$padj, breaks=100))
dev.off()

resSig <- res[res$padj <0.05,]
head (resSig[order(resSig$pval),])
write.csv(resSig, file="genome_differential_expression_EZH1_Mut_single.csv")
write.csv(res, file="genome_assayed_genes_EZH1_Mut_single.csv")

cdsBlind <- estimateDispersions (cds, method="blind")
vsd <- getVarianceStabilizedData (cdsBlind)
select <-order(res$padj)[1:50]
colors<- colorRampPalette(c("white", "darkblue"))(100)
pdf("genome_heatmap_EZH1_Mut_single.pdf")
print(heatmap(vsd[select,], col=colors, scale="none"))
dev.off()

res = nbinomTest (cds, "WT", "ezh2_ezh1")

pdf("plotMA_EZH2_EZH1_Mut_single.pdf")
print(plotMA(res))
dev.off()

pdf("histogram_EZH2_EZH1_Mut_single.pdf")
print(hist(res$padj, breaks=100))
dev.off()

resSig <- res[res$padj <0.05,]
head (resSig[order(resSig$pval),])
write.csv(resSig, file="genome_differential_expression_EZH2_EZH1_Mut_single.csv")
write.csv(res, file="genome_assayed_genes_EZH2_EZH1_Mut_single.csv")

cdsBlind <- estimateDispersions (cds, method="blind")
vsd <- getVarianceStabilizedData (cdsBlind)
select <-order(res$padj)[1:50]
colors<- colorRampPalette(c("white", "darkblue"))(100)
pdf("genome_heatmap_EZH2_EZH1_Mut_single.pdf")
print(heatmap(vsd[select,], col=colors, scale="none"))
dev.off()

res = nbinomTest (cds, "WT", "eed_het")

pdf("plotMA_EED_HET_Mut_single.pdf")
print(plotMA(res))
dev.off()

pdf("histogram_EED_HET_Mut_single.pdf")
print(hist(res$padj, breaks=100))
dev.off()

resSig <- res[res$padj <0.05,]
head (resSig[order(resSig$pval),])
write.csv(resSig, file="genome_differential_expression_EED_HET_Mut_single.csv")
write.csv(res, file="genome_assayed_genes_EED_HET_Mut_single.csv")

cdsBlind <- estimateDispersions (cds, method="blind")
vsd <- getVarianceStabilizedData (cdsBlind)
select <-order(res$padj)[1:50]
colors<- colorRampPalette(c("white", "darkblue"))(100)
pdf("genome_heatmap_EED_HET_Mut_single.pdf")
print(heatmap(vsd[select,], col=colors, scale="none"))
dev.off()

res = nbinomTest (cds, "WT", "eed")

pdf("plotMA_EED_Mut_single.pdf")
print(plotMA(res))
dev.off()

pdf("histogram_EED_Mut_single.pdf")
print(hist(res$padj, breaks=100))
dev.off()


resSig <- res[res$padj <0.05,]
head (resSig[order(resSig$pval),])
write.csv(resSig, file="genome_differential_expression_EED_Mut_single.csv")
write.csv(res, file="genome_assayed_genes_EED_Mut_single.csv")

cdsBlind <- estimateDispersions (cds, method="blind")
vsd <- getVarianceStabilizedData (cdsBlind)
select <-order(res$padj)[1:50]
colors<- colorRampPalette(c("white", "darkblue"))(100)
pdf("genome_heatmap_EED_Mut_single.pdf")
print(heatmap(vsd[select,], col=colors, scale="none"))
dev.off()
