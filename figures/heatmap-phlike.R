rm(list=ls())

ann <- read.delim("/mnt/projects/ikaros/data/samples.csv", stringsAsFactors = F)
rownames(ann) <- ann$UPN
ann <- ann[ann$Xeno != "yes" & ann$Exclude != "yes",]

# same order of samples as in figure 4A of paper
ann <- ann[c("92D", "N7R", "GI8R", "DL2R", "HV57R", "GI13R", "HV80D", "HV80R", "841D", "400D", "1060D", "802D", "887D", "506D", "961D", "BB16D", "839R", "1089D", "GI8D", "AL9890R", "360D", "365D"),]
ann$IKZF1 <- ann$Group
ann$IKZF1[ann$IKZF1=="IK6"] <- "IKN"
ann$IKZF1[ann$Alias=="IKD_5"] <- "IKN"
ann$IKZF1 <- factor(ann$IKZF1, levels=c("IKN", "IKD", "IKC"))
ann$RASpw <- factor(ifelse(ann$RAS == "-", "wt", "mut"))
ann$JAKpw <- factor(ifelse(ann$JAK_STAT == "-", "wt", "mut"))

annotation.colors <- list(
  IKZF1=c("IKN" = "red", "IKD" = "blue", "IKC" = "gray"),
  PC=c("yes"="red", "no"="blue"),
  Relapsing=c("yes"="red", "no"="blue"),
  Timepoint=c("diagnosis" = "blue", "relapse" = "red"),
  RASpw=c("wt"="blue", "mut"="red"),
  JAKpw=c("wt"="blue", "mut"="red")
)

m <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqExprMatrix/table.csv", stringsAsFactors = F, row.names = 1)
m <- as.matrix(m)
m <- m[,ann$Alias]
colnames(m) <- ann$UPN[match(colnames(m), ann$Alias)]

# Ph-like signature from Roberts et al. (2012), Table S2 (http://www.sciencedirect.com/science/article/pii/S1535610812002541)
library("biomaRt")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75
hgnc <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), mart=mart)
hgnc <- hgnc[hgnc$hgnc_symbol != "" & grepl("ENSG", hgnc$ensembl_gene_id),]
mullighan.tableS15 <- read.delim("/mnt/projects/ikaros/data/public/harvey_mullighan_2010_suppTable15_Top 100 Rank Order Genes Defining ROSE Cluster R8.txt", stringsAsFactors = F)
mullighan.tableS15.ensembl <- merge(data.frame(hgnc=unique(mullighan.tableS15$Symbol[!is.na(mullighan.tableS15$Symbol) & mullighan.tableS15$Symbol != ""])), hgnc, by.x="hgnc", by.y="hgnc_symbol")

roberts.2012 <- read.delim("/mnt/projects/ikaros/data/public/roberts_2012_tableS2.Ph-like.vs.other.txt", stringsAsFactors = F, check.names = F)
roberts.2012 <- roberts.2012[roberts.2012$Gene_Symbol != "---", c("Gene_Symbol", "FoldChange(G2/G1)", "FDR")]
names(roberts.2012) <- c("hgnc", "log2FC", "q")
roberts.2012 <- roberts.2012[order(roberts.2012$q),]
roberts.2012 <- roberts.2012[!duplicated(roberts.2012$hgnc),]
roberts.2012 <- merge(roberts.2012, hgnc, by.x="hgnc", by.y="hgnc_symbol")
roberts.2012 <- roberts.2012[roberts.2012$ensembl_gene_id %in% rownames(m),]
roberts.2012 <- roberts.2012[order(roberts.2012$q),]
roberts.2012.top50 <- rbind(roberts.2012[roberts.2012$log2FC>0,][1:50,],roberts.2012[roberts.2012$log2FC<0,][1:50,])
roberts.2012.top50$hgnc[roberts.2012.top50$log2FC<0] <- paste0(roberts.2012.top50$hgnc[roberts.2012.top50$log2FC<0], "*")

gene.sets <- list(
  "mullighan.patent" = read.delim("/mnt/projects/ikaros/results/qlucore/mullighan_patent_phlike-15.qlucore.txt"),
  "mullighan.2010.tableS15" = read.delim("/mnt/projects/ikaros/results/qlucore/harvey_mullighan_2010_suppTable15_Top 100 Rank Order Genes Defining ROSE Cluster R8.qlucore.txt", stringsAsFactors = F),
  "roberts.2012.tableS2.top50.UPandDN" = roberts.2012.top50
)
# remove duplicate entry for LAIR1
gene.sets[["mullighan.2010.tableS15"]] <- gene.sets[["mullighan.2010.tableS15"]][gene.sets[["mullighan.2010.tableS15"]]$ensembl_gene_id != "ENSG00000262936",]

pdf("/mnt/projects/p2ry8-crlf2/results/figures/expression-heatmap-phlike-samples-unclustered.pdf", height=12)
for (setname in names(gene.sets)) {
  mf <- m[rownames(m) %in% gene.sets[[setname]]$ensembl_gene_id,]
  rownames(mf) <- gene.sets[[setname]]$hgnc
  dist.genes <- as.dist(1-cor(t(mf), method = "spearman"))
  dist.samples <- dist(t(mf), method = "euclidean")
  
  library(pheatmap)
  pheatmap(
    mf, 
    main = setname,
    cluster_cols = F,
    clustering_method = "average",
    clustering_distance_rows = dist.genes, 
    clustering_distance_cols = dist.samples, 
    annotation_col = ann[,c("IKZF1", "PC", "Relapsing", "Timepoint", "RASpw", "JAKpw"),drop=F],
    scale="row",
    cellheight = 7,
    legend = F,
    fontsize = 7,
    #  labels_col = ann$UPN[match(attr(dist.samples, "Labels"), rownames(sample.table))],
    col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
    annotation_colors = annotation.colors
  )
}
dev.off()

