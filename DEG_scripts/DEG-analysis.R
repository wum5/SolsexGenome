# load relevant R packages
library(edgeR)
library(statmod)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(rgl)

# read data
# 3 male individuals: 670-1190, 670-105, 716-6 (670 and 716 are two populations)
# 3 female individuals: 670-34, 716-23, 716-4 
# tissues: B represents buds while F represents flowers
rawdata <- read.delim("raw_featurecount_data.txt", header = TRUE)
Pop <- factor(c("1","1","1","1","1","1","2","2","2","2","2","2"))
Tissue <- factor(c("B","F","B","F","B","F","B","F","B","F","B","F"))
Sex <- factor(c("M","M","M","M","F","F","F","F","F","F","M","M"))

# filtering and normalization
y <- DGEList(counts=rawdata[,3:ncol(rawdata)], genes=rawdata[,1:2])
d <- duplicated(y$genes$Geneid)
y <- y[!d,]
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples

# multidimensional scaling plots
# here we can see the gene expression are firstly separated by tissues, then by populations
# in this case, sex is a minor factor for gene differential expression
pdf("MDSplot_rna.pdf", width=6, height=6)
col.cell <- c("red", "orange", "blue", "cyan")
mds <- plotMDS(y$counts, col=col.cell, pch=c(19,15)[Sex], ndim=2)
legend("topright", fill = c("red", "orange", "blue", "cyan"), 
       legend = c("670-Bud", "670-flower", "716-Bud", "716-flower"))
dev.off()

# variances being explained by 
s <- svd(log(y$counts+1,2)-log(rowMeans(y$counts+1,2)))
plot(s$d^2/sum(s$d^2))



# experimental desgin (here we used a GLM considering the effect of population, tissue, 
# sex and the interaction between tissue and sex)
design <- model.matrix(~Pop+Tissue+Sex+Tissue:Sex)
rownames(design) <- colnames(y)

# estimating the dispersion
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)

# differential expression
fit <- glmFit(y, design)
lrt <- glmLRT(fit,4)
plotMD(lrt, status=decideTestsDGE(lrt))
abline(h=c(-1, 1), col="blue")

# number of DEGs (by using the GLM, not too many DE genes were identified)
summary(decideTestsDGE(lrt))
DEG <- summary(decideTestsDGE(lrt))
DEG_num <- (DEG[1]+DEG[3])

# heatmap of the DEG 
# here the heatmap shows our model can pull out sex-biased genes (they are separated)
o <- order(lrt$table$PValue)
DEG_list <- cpm(y)[o[1:DEG_num],]
geneValsLog = log(DEG_list + 1)
pheatmap(geneValsLog, treeheight_row = 0, show_rownames = F)

# generate ouput files and plots
top100DEG <- topTags(lrt, n=100)
DEgenes <- top100DEG$table[top100DEG$table$FDR<0.05 & abs(top100DEG$table$logFC)>1,]
write.csv(DEgenes, file="DEgenes_by_GLM.csv", row.names = F)


# To indentify more DE genes, I used a pairwise comparison (i.e. male vs. femlae) 
# to pull out genes differential expressed in at least one of the two tissues
# First in flower buds
# filtering and normalization 
y <- DGEList(counts=rawdata[,c(3,5,7,9,11,13)], genes=rawdata[,1:2])
d <- duplicated(y$genes$Geneid)
y <- y[!d,]
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples

Pop <- factor(c("1","1","1","2","2","2"))
Sex <- factor(c("M","M","F","F","F","M"))

design <- model.matrix(~Sex)
rownames(design) <- colnames(y)

# estimating the dispersion
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)

# differential expression
fit <- glmFit(y, design)
lrt <- glmLRT(fit,2)

# very few genes were detected to be differentially expressed in flower buds
plotMD(lrt, status=(decideTestsDGE(lrt)), main="flower buds", 
       hl.col=c("green3", "red"),ylim=c(-9.2,9.2))
abline(h=c(-1, 1), col="blue")

# generate ouput files and plots
topTags(lrt)
summary(decideTestsDGE(lrt))

top100DEG <- topTags(lrt,100)
write.csv(topTags(lrt, 100000), file="allGene_Expr_bud.csv", row.names = F)

bud_DEG <- top100DEG$table[top100DEG$table$FDR<0.05 & abs(top100DEG$table$logFC)>1,]
write.csv(bud_DEG, file="bud_DEGs.csv", row.names = F)


# Second in mature flower 
# filtering and normalization 
y <- DGEList(counts=rawdata[,c(4,6,8,10,12,14)], genes=rawdata[,1:2])
d <- duplicated(y$genes$Geneid)
y <- y[!d,]
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples

Pop <- factor(c("1","1","1","2","2","2"))
Sex <- factor(c("M","M","F","F","F","M"))

design <- model.matrix(~Sex)
rownames(design) <- colnames(y)

# estimating the dispersion
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)

# differential expression
fit <- glmFit(y, design)
lrt <- glmLRT(fit,2)

# ~100 genes were detected to be differentially expressed in mature flowers
plotMD(lrt, status=(decideTestsDGE(lrt)),  main="mature flower", ylim=c(-9.2,9.2))
abline(h=c(-1, 1), col="blue")

# generate ouput files and plots
topTags(lrt)
summary(decideTestsDGE(lrt))

top100DEG <- topTags(lrt,100)
write.csv(topTags(lrt, 100000), file="allGene_Expr_flower.csv", row.names = F)

flower_DEG <- top100DEG$table[top100DEG$table$FDR<0.05 & abs(top100DEG$table$logFC)>1,]
write.csv(flower_DEG, file="flower_DEGs.csv", row.names = F)


# check how similarity between the DEGs in the two tissues
flower_DEGid <- flower_DEG$Geneid
bud_DEGid <- bud_DEG$Geneid
shared_DEGid <- intersect(flower_DEGid, bud_DEGid) 

# number of DEGs in mature flowers
length(flower_DEGid)

# number of DEGs in flower buds
length(bud_DEGid)

# number of DEGs in both tissues
length(shared_DEGid)
shared_DEGid


# heatmap of the DEGs (combining DEGs in the two tissues)
# here the heatmap shows our model can pull out sex-biased genes (they are separated)
y <- DGEList(counts=rawdata[,3:ncol(rawdata)], genes=rawdata[,1:2])
d <- duplicated(y$genes$Geneid)
y <- y[!d,]
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

combined_DEGid <- union(flower_DEGid, bud_DEGid)
DEG_counts <- y$counts[y$genes$Geneid %in% combined_DEGid, ]
geneValsLog = log(DEG_counts + 1)
pheatmap(geneValsLog, treeheight_row = 0, show_rownames = F)

write.csv(combined_DEGid, file="DEG_by_pairwiseCompr.csv", row.names = F)

