# read sample table

sample.table <- read.delim("/mnt/projects/ikaros/data/all_samples_unsupervised_clustering.csv", stringsAsFactors = F, row.names = "ID", na.strings = "")
sample.table$IKZF1 <- as.factor(sample.table$IKZF1)
sample.table$Subtype[grepl("CD19", sample.table$Subtype)] <- "CD19"

# heatmap annotations 

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

# biomart for gene names

library(biomaRt)
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl") # GRCh37, v75
hgnc <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), mart=mart)

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
counts <- counts[rowSums(counts) >= 100,] ; dim(counts)

# EDAseq within-lane (GC content) and between-lane (library size) normalization

gclength <- getBM(attributes=c("ensembl_gene_id", "percentage_gc_content", "transcript_length"), mart=mart)
gclength <- gclength[order(gclength$ensembl_gene_id, -gclength$transcript_length),]
gclength <- gclength[!duplicated(gclength$ensembl_gene_id),]
gclength <- data.frame(id=gclength$ensembl_gene_id, gc=gclength$percentage_gc_content, length=gclength$transcript_length)
rownames(gclength) <- gclength$id
gclength <- gclength[,c("gc", "length")]

library(EDASeq)
common <- intersect(rownames(counts), rownames(gclength))
counts <- counts[common,] ; dim(counts)
gclength <- gclength[common,]
edaseq.counts <- newSeqExpressionSet(counts=counts, featureData = gclength)
edaseq.withinLane <- withinLaneNormalization(edaseq.counts, "gc", which="full", offset=TRUE, round=FALSE)
edaseq.betweenLane <- betweenLaneNormalization(edaseq.withinLane, which="full", offset=TRUE, round=FALSE)

# replace ensembl IDs with HGNC gene names in count matrix

rownames(counts) <- ifelse(hgnc$hgnc_symbol[match(rownames(counts), hgnc$ensembl_gene_id)] != "", hgnc$hgnc_symbol[match(rownames(counts), hgnc$ensembl_gene_id)], rownames(counts))

# DESeq2

subtypes <- sample.table[, "Subtype", drop=F]
subtypes$Subtype[sample.table$UPN %in% c("N7R", "AL9890R", "GI8D", "HV57D", "GI8R", "HV57R", "DL2R")] <- "excluded"
subtypes$Subtype <- as.factor(subtypes$Subtype)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = subtypes, design = ~Subtype)
normalizationFactors(dds) <- exp(-1 * offst(edaseq.betweenLane))
dds <- DESeq(dds, minReplicatesForReplace=5)
vst <- assay(varianceStabilizingTransformation(dds, blind=TRUE))

# write out gene expression matrix
write.table(
  vst, 
  "/mnt/projects/p2ry8-crlf2/results/ALL-primary-normalizedExpressionMatrix-for-unsupervised-clustering.tsv", 
  col.names = paste0(sample.table$UPN[match(colnames(vst), rownames(sample.table))], "(", sample.table$Subtype[match(colnames(vst), rownames(sample.table))], ")"), 
  row.names = T, sep="\t", quote = F
)

# get DEGs for all pairwise comparisons

results <- list()
for (comparison in c("iAMP21.vs.PC-only", 
                     "iAMP21.vs.ER", 
                     "iAMP21.vs.HD", 
                     "iAMP21.vs.B-other", 
                     "PC-only.vs.ER", 
                     "PC-only.vs.HD", 
                     "PC-only.vs.B-other", 
                     "ER.vs.HD", 
                     "ER.vs.B-other", 
                     "HD.vs.B-other",
                     "CD19.vs.iAMP21",
                     "CD19.vs.PC-only",
                     "CD19.vs.ER",
                     "CD19.vs.HD")) {
  groups <- unlist(strsplit(comparison, ".vs."))
  groupA <- groups[1]
  groupB <- groups[2]
  print(paste("Comparing", groupA, "with", groupB, "..."))
  res <- as.data.frame(results(dds, contrast=c("Subtype", groupA, groupB), cooksCutoff=FALSE))
  res <- res[order(res$padj),]
  results[[comparison]] <- res
}

# gene filtering
degs <- character(0)
topN <- 10
maxQ <- 0.001
minFC <- 0.5
for (comparison in names(results)){
  res <- results[[comparison]]
  
  sigUp <- res[!is.na(res$padj) & res$padj <= maxQ & !is.na(res$log2FoldChange) & res$log2FoldChange >= minFC,]
  sigUp <- sigUp[!grepl("^ENSG0", rownames(sigUp)),]
  sigUp <- sigUp[1:min(nrow(sigUp), topN),]

  sigDn <- res[!is.na(res$padj) & res$padj <= maxQ & !is.na(res$log2FoldChange) & res$log2FoldChange <= -minFC,]
  sigDn <- sigDn[!grepl("^ENSG0", rownames(sigDn)),]
  sigDn <- sigDn[1:min(nrow(sigDn), topN),]
  
  degs <- c(degs, rownames(sigUp))
  degs <- c(degs, rownames(sigDn))
}
degs <- unique(degs)

# heatmap

mf <- vst[degs,]
mf <- t(scale(t(mf)))
mf[mf<=-4] <- -4; mf[mf>=4] <- 4

pdf(paste0("/mnt/projects/p2ry8-crlf2/results/rnaseq-clustering-top", topN, "Genes-pairwise-comparisons.pdf"), width=10, height=10)

library(pheatmap)
pheatmap(
  mf,
  main = paste("Top", topN, "up/dn-regulated genes from each pairwise comparisons of subtypes"),
  clustering_method = "average",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "euclidean",
  col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
  fontsize = 7,
  fontsize_row = 3,
  border_color=NA,
  labels_col = sample.table$UPN[match(colnames(mf), rownames(sample.table))],
  annotation_col = sample.table[,c("Subtype", "IKZF1", "PC", "chr21", "RNAExtrDate", "RNAKit", "PrepDate", "Flowcell", "Origin", "Timepoint"),drop=F],
  annotation_colors = annotation.colors,
  legend=F
)

dev.off()
