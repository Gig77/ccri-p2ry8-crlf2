# read sample annotation

sample.table <- read.delim("/mnt/projects/ikaros/data/all_samples_unsupervised_clustering.csv", stringsAsFactors = F, row.names = "Alias", na.strings = "")
sample.table <- sample.table[!sample.table$Subtype %in% c("MLL pos", "BCR-ABL1"),]
sample.table <- sample.table[!grepl("CD19", sample.table$Subtype),]
sample.table <- sample.table[!(sample.table$Subtype %in% c("HD") & sample.table$chr21 %in% c("DS")),]
sample.table <- sample.table[!sample.table$Subtype %in% c("PC-only") | sample.table$chr21 %in% c("DS"),]
sample.table$Subtype[sample.table$Subtype=="PC-only"] <- "DS"
sample.table$Subtype <- factor(sample.table$Subtype)

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
counts <- counts[rowSums(counts) >= 10,]

boxplot(log2(counts+0.1))

# EDAseq within-lane and between-lane normalization

library(biomaRt)
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl") # GRCh37, v75
hgnc <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), mart=mart)

gclength <- getBM(attributes=c("ensembl_gene_id", "percentage_gc_content", "transcript_length"), mart=mart)
gclength <- gclength[order(gclength$ensembl_gene_id, -gclength$transcript_length),]
gclength <- gclength[!duplicated(gclength$ensembl_gene_id),]
gclength <- data.frame(id=gclength$ensembl_gene_id, gc=gclength$percentage_gc_content, length=gclength$transcript_length)
rownames(gclength) <- gclength$id
gclength <- gclength[,c("gc", "length")]

library(EDASeq)
common <- intersect(rownames(counts), rownames(gclength))
edaseq.counts <- newSeqExpressionSet(counts=counts[common,], featureData = gclength[common,])
edaseq.withinLane <- withinLaneNormalization(edaseq.counts, "gc", which="full", offset=TRUE, round=FALSE)
edaseq.betweenLane <- betweenLaneNormalization(edaseq.withinLane, which="full", offset=TRUE, round=FALSE)

#counts.norm <- normCounts(edaseq.betweenLane)
#mode(counts.norm) <- "integer"
#boxplot(log2(counts.norm+0.1))

# DESeq2

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts[common,], colData = sample.table[, "Subtype", drop=F], design = ~Subtype)
normalizationFactors(dds) <- exp(-1 * offst(edaseq.betweenLane))
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
dds <- replaceOutliers(dds, trim = 0.2, cooksCutoff = 0.7, minReplicates = 7)
dds <- nbinomWaldTest(dds)

save(dds, file="/mnt/projects/p2ry8-crlf2/results/ALL-primary-forFileMaker.DESeq2.Rdata")

res.hd <- as.data.frame(results(dds, contrast=c("Subtype", "ER", "HD"), cooksCutoff=FALSE))
res.hd <- res.hd[order(res.hd$padj),]

res.iamp <- as.data.frame(results(dds, contrast=c("Subtype", "ER", "iAMP21"), cooksCutoff=FALSE))
res.iamp <- res.iamp[order(res.iamp$padj),]

res.ds <- as.data.frame(results(dds, contrast=c("Subtype", "ER", "DS"), cooksCutoff=FALSE))
res.ds <- res.ds[order(res.ds$padj),]

res.bother <- as.data.frame(results(dds, contrast=c("Subtype", "ER", "B-other"), cooksCutoff=FALSE))
res.bother <- res.bother[order(res.bother$padj),]
