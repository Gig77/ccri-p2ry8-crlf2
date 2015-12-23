library(ExomeDepth)

exons.bed <- "/mnt/projects/generic/data/illumina/nexterarapidcapture_exome_targetedregions.nochr.bed"
bam.files <- list.files(path="/mnt/projects/p2ry8-crlf2/data/bam", pattern="*.bam$", full.names=T)
bam.files <- bam.files[!grepl("(abra|715C|715D|715C|715R_)", bam.files)]
reference.file <- "/mnt/projects/generic/data/broad/human_g1k_v37.fasta"

# get counts
counts <- getBamCounts(bed.file = exons.bed, bam.files = bam.files, include.chr = FALSE, referenceFasta = reference.file)

# rename columns
sample.names <- paste0(sub(".*variant_calling_process_sample_(.+)_realigned.*", "\\1", bam.files))
names(counts@values@unlistData)[-1] <- sample.names

# save
save(counts, file="/mnt/projects/p2ry8-crlf2/results/exomeDepth/counts.RData")
