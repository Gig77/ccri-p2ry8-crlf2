rm(list=ls())

ann <- read.delim("/mnt/projects/ikaros/data/samples.csv", stringsAsFactors = F)
ann <- ann[ann$Xeno != "yes" & ann$Exclude != "yes",]
ann$Group[ann$Group=="IK6"] <- "IKN"
ann$Group[ann$Alias=="IKD_5"] <- "IKN"
ann$Group <- factor(ann$Group, levels=c("IKN", "IKD", "IKC"))

m <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqExprMatrix/table.csv", stringsAsFactors = F, row.names = 1)
m <- as.matrix(m)
m <- m[,colnames(m) %in% ann$Alias]
colnames(m) <- ann$UPN[match(colnames(m), ann$Alias)]

topN <- 50

IKNp.vs.IKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKNp_vs_IKCp/table.csv", stringsAsFactors = F)
IKNp.vs.IKCp.filt <- IKNp.vs.IKCp[!is.na(IKNp.vs.IKCp$p) & IKNp.vs.IKCp$p <= 1e-8,]
IKNp.vs.IKCp.filt <- IKNp.vs.IKCp.filt[order(IKNp.vs.IKCp.filt$fc, decreasing=T),]
IKNp.vs.IKCp.filt <- IKNp.vs.IKCp.filt[c(1:topN,(nrow(IKNp.vs.IKCp.filt)-topN+1):nrow(IKNp.vs.IKCp.filt)),]

IKD.vs.IKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKD_vs_IKCp/table.csv", stringsAsFactors = F)
IKD.vs.IKCp.filt <- IKD.vs.IKCp[!is.na(IKD.vs.IKCp$p) & IKD.vs.IKCp$p <= 2e-3,]
IKD.vs.IKCp.filt <- IKD.vs.IKCp.filt[order(IKD.vs.IKCp.filt$fc, decreasing=T),]
IKD.vs.IKCp.filt <- IKD.vs.IKCp.filt[unique(c(1:topN,(nrow(IKD.vs.IKCp.filt)-topN+1):nrow(IKD.vs.IKCp.filt))),]

# filter expression matrix for differentially regulated genes
mf <- m[rownames(m) %in% c(IKNp.vs.IKCp.filt$ids, IKD.vs.IKCp.filt$ids),]
hgnc <- IKNp.vs.IKCp$Gene[match(rownames(mf), IKNp.vs.IKCp$ids)]
rownames(mf) <- ifelse(is.na(hgnc), rownames(mf), hgnc)
hc.genes <- hclust(as.dist(1-cor(t(mf), method = "spearman")), method = "average")
#hc.samples <- hclust(dist(1-cor(mf)), method = "average")
hc.samples <- hclust(dist(t(mf), method = "euclidean"))
plot(hc.samples)

library(ggplot2)
pdf("/mnt/projects/p2ry8-crlf2/results/figures/expression-heatmap.pdf")
heatmap.2(
  mf, 
  main="Heatmap", 
  Rowv=as.dendrogram(hc.genes),
  Colv=reorder(as.dendrogram(hc.samples), ann$Group),
  colsep=seq(1:dim(mf)[[2]]), 
  #rowsep=seq(1:dim(mf)[[1]]), 
  sepcolor="grey92", 
  sepwidth=c(0.005,0.005),
  ColSideColors=c("red", "blue", "gray")[ann$Group], 
  scale="row",
  keysize=0.8, 
  col=rev(redblue(256)), 
  density.info="none", 
  symkey=FALSE, 
  labRow=rownames(mf), 
  labCol=colnames(mf), 
  trace="none",
  cex.main=0.7, 
  cexRow=0.3, 
  cexCol=1, 
  mar=c(5,5) 
)
dev.off()
