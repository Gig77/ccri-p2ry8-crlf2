library(biomaRt)
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl") # GRCh37, v75
hgnc <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), mart=mart)

# read sample annotation

sample.table <- read.delim("/mnt/projects/ikaros/data/all_samples_unsupervised_clustering.csv", stringsAsFactors = F, row.names = "ID", na.strings = "")
sample.table$IKZF1 <- as.factor(sample.table$IKZF1)
sample.table$Subtype[grepl("CD19", sample.table$Subtype)] <- "CD19"

# build count matrix

counts <- NULL
for (i in 1:nrow(sample.table)) {
  sample <- sample.table[i,]
  if (is.null(counts)) {
    counts <- read.delim(sample$File, row.names = 1, header = F, comment.char = "_")
    names(counts) <- rownames(sample)
  } else {
    counts.sample <- read.delim(sample$File, row.names = 1, header = F, comment.char = "_")
    counts[,rownames(sample)] <- counts.sample[rownames(counts),]
  }
}
counts <- as.matrix(counts)
counts[is.na(counts)] <- 0

# filter lowly expressed genes

counts.filtered <- counts[rowSums(counts) >= 100,]
print(sprintf("Clustering based on %d expressed genes", nrow(counts.filtered)))

# EDAseq within-lane and between-lane normalization

gclength <- getBM(attributes=c("ensembl_gene_id", "percentage_gc_content", "transcript_length"), mart=mart)
gclength <- gclength[order(gclength$ensembl_gene_id, -gclength$transcript_length),]
gclength <- gclength[!duplicated(gclength$ensembl_gene_id),]
gclength <- data.frame(id=gclength$ensembl_gene_id, gc=gclength$percentage_gc_content, length=gclength$transcript_length)
rownames(gclength) <- gclength$id
gclength <- gclength[,c("gc", "length")]

library(EDASeq)
common <- intersect(rownames(counts.filtered), rownames(gclength))
edaseq.counts <- newSeqExpressionSet(counts=counts.filtered[common,], featureData = gclength[common,])
edaseq.withinLane <- withinLaneNormalization(edaseq.counts, "gc", which="full", offset=TRUE, round=FALSE)
edaseq.betweenLane <- betweenLaneNormalization(edaseq.withinLane, which="full", offset=TRUE, round=FALSE)
#counts.edaseq.log <- log(normCounts(edaseq.betweenLane) + 0.1, 2)

# DESeq2 variant-stabilizing transformation using EDASeq normalization factors

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts.filtered[common,], colData = sample.table[, "Subtype", drop=F], design = ~1)
normalizationFactors(dds) <- exp(-1 * offst(edaseq.betweenLane))
vst <- varianceStabilizingTransformation(dds, blind=TRUE)
counts.normalized.transformed <- assay(vst)

pdf("/mnt/projects/p2ry8-crlf2/results/rnaseq-unsupervised-clustering.pdf", width=10, height=10)

# ----
# sample-by-sample heatmap
# ----

# filter for most variable genes across samples

percentMostVariable <- 0.4
cov <- apply(counts.normalized.transformed, 1, function(x) (quantile(x, 0.75)-quantile(x, 0.25))/quantile(x, 0.5))
counts.normalized.transformed.mostvariable <- counts.normalized.transformed[cov >= quantile(cov, 1-percentMostVariable),]

# get distance matrix, mask diagonal for improved color range

dists <- as.dist(1-cor(counts.normalized.transformed.mostvariable, method="spearman")) # spearman correlation
dists.matr <- as.matrix(dists)
diag(dists.matr) <- rowMeans(dists.matr)

annotation.colors <- list(
  Subtype=c("CD19" = "gray", "ER" = "red", "HD" = "blue", "iAMP21" = "yellow", "PC-only" = "black", "B-other" = "lightblue", "MLL pos" = "lightgreen", "BCR-ABL1" = "brown"),
  IKZF1=c("del" = "blue", "dn" = "red", "wt" = "gray"),
  PC=c("yes" = "black", "no" = "gray"),
  chr21=c("normal" = "gray", "DS" = "blue", "trisomy" = "red", "iAMP21" = "yellow"),
  RNAExtrDate=c("2014-02-10" = "gray", "2014-08-01" = "red", "2014-09-30" = "blue", "2015-04-02" = "yellow", "2015-04-07" = "black", "2015-08-28" = "lightblue", "Kamilla" = "lightgreen", "2016-02-03" = "orange", "2016-02-09" = "brown"),
  RNAKit=c("Zymo" = "gray", "mirVana" = "black"),
  PrepDate=c("2014-09-16" = "gray", "2014-09-18" = "red", "2014-10-13" = "blue", "2015-05-05" = "yellow", "2015-05-06" = "black", "2016-02-09/10" = "lightblue", "2016-02-11/12" = "lightgreen", "2016-02-05/08" = "orange"),
  Flowcell=c("C57C3ACXX" = "gray", "C7K1TACXX" = "blue", "C88A8ACXX" = "red", "C97W1ANXX" = "yellow"),
  Origin=c("AT" = "gray", "GER" = "black"),
  Timepoint=c("diagnosis" = "gray", "relapse" = "black")
)

library(pheatmap)
pheatmap(
  dists.matr, 
  main = sprintf("Clustering based on %d genes\n(%d%% most variable expressed genes)", nrow(counts.normalized.transformed.mostvariable), as.integer(percentMostVariable*100)),
  clustering_method = "average",
  clustering_distance_rows = dists, 
  clustering_distance_cols = dists, 
  annotation_col = sample.table[,c("Subtype", "IKZF1", "PC", "chr21", "RNAExtrDate", "RNAKit", "PrepDate", "Flowcell", "Origin", "Timepoint"),drop=F],
  legend = F,
  fontsize = 6,
  labels_col = sample.table$UPN[match(attr(dists, "Labels"), rownames(sample.table))],
  col = colorRampPalette(brewer.pal(10, "RdBu"))(256),
  annotation_colors = annotation.colors
)

# ----
# gene-by-sample heatmaps
# ----

genesets <- list()

# most differentially expressed genes from IKN comparison

topN <- 30
IKNp.vs.IKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKNp_vs_IKCp/table.csv", stringsAsFactors = F)
IKNp.vs.IKCp.filt <- IKNp.vs.IKCp[!is.na(IKNp.vs.IKCp$p) & IKNp.vs.IKCp$p <= 1e-8,]
IKNp.vs.IKCp.filt <- IKNp.vs.IKCp.filt[order(IKNp.vs.IKCp.filt$fc, decreasing=T),]
IKNp.vs.IKCp.filt <- IKNp.vs.IKCp.filt[c(1:topN,(nrow(IKNp.vs.IKCp.filt)-topN+1):nrow(IKNp.vs.IKCp.filt)),]
genesets[[sprintf("Top-%d DEGs IKNp-vs-IKCp", topN*2)]] <- IKNp.vs.IKCp.filt$ids

# MSigDB gene sets

gmt <- read.delim("/mnt/projects/generic/data/msigdb5.0/msigdb.v5.0.symbols.gmt", stringsAsFactors = F, header=F)
for (geneset in c("PID_FAK_PATHWAY", "PID_VEGFR1_2_PATHWAY", "KEGG_JAK_STAT_SIGNALING_PATHWAY")) {
  geneset.genes <- gmt[toupper(gmt$V1)==geneset,c(-1,-2)]
  geneset.genes <- sort(unique(geneset.genes[geneset.genes != ""]))
  geneset.genes <- hgnc$ensembl_gene_id[match(geneset.genes, hgnc$hgnc_symbol)]
  geneset.genes <- geneset.genes[!is.na(geneset.genes)]
  genesets[[geneset]] <- geneset.genes
}

# DSigDB gene sets

gmt <- read.delim("/mnt/projects/generic/data/DSigDB/DSigDB_v1.0_All.gmt", stringsAsFactors = F, header=F)
for (geneset in c("DASATINIB_FDA")) {
  geneset.genes <- gmt[toupper(gmt$V1)==geneset,c(-1,-2)]
  geneset.genes <- sort(unique(geneset.genes[geneset.genes != ""]))
  geneset.genes <- hgnc$ensembl_gene_id[match(geneset.genes, hgnc$hgnc_symbol)]
  geneset.genes <- geneset.genes[!is.na(geneset.genes)]
  genesets[[geneset]] <- geneset.genes
}

# curated gene sets
gmt <- read.delim("/mnt/projects/ikaros/data/ikaros_curated_genesets_gsea.gmt", stringsAsFactors = F, header=F)
for (geneset in gmt$V1) {
  geneset.genes <- gmt[gmt$V1==geneset,c(-1,-2)]
  geneset.genes <- sort(unique(geneset.genes[geneset.genes != ""]))
  geneset.genes <- hgnc$ensembl_gene_id[match(geneset.genes, hgnc$hgnc_symbol)]
  geneset.genes <- geneset.genes[!is.na(geneset.genes)]
  genesets[[geneset]] <- geneset.genes
}

genesets[["IACOBUCCI_2012_IKAROS_DELETED_UP_OR_DN"]] <- c(genesets[["IACOBUCCI_2012_IKAROS_DELETED_UP"]], genesets[["IACOBUCCI_2012_IKAROS_DELETED_DN"]])

# plot clustered heatmap for each gene set

for (geneset in names(genesets)) {
  print(geneset)
  geneset.genes <- genesets[[geneset]]
  mf <- counts.edaseq.log[rownames(counts.edaseq.log) %in% geneset.genes,]
  mf <- mf[apply(mf, 1, sd) > 0,]
  rownames(mf) <- ifelse(hgnc$hgnc_symbol[match(rownames(mf), hgnc$ensembl_gene_id)] != "", hgnc$hgnc_symbol[match(rownames(mf), hgnc$ensembl_gene_id)], rownames(mf))

  pheatmap(
    mf,
    main = geneset,
    clustering_method = "average",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "euclidean",
    col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
    scale="row",
    fontsize = 6,
    fontsize_row = 5,
    labels_col = sample.table$UPN[match(colnames(mf), rownames(sample.table))],
    annotation_col = sample.table[,c("Subtype", "IKZF1", "PC", "chr21", "RNAExtrDate", "RNAKit", "PrepDate", "Flowcell", "Origin", "Timepoint"),drop=F],
    annotation_colors = annotation.colors,
    legend=F
  )
}

dev.off()
