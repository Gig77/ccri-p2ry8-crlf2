library(GenomicRanges)
library(ExomeDepth)
library(optparse)

# read/check parameters
option_list <- list(make_option("--sample", type="character", help="sample name"))
opt <- parse_args(OptionParser(option_list=option_list))
#opt <- data.frame(sample="HV80D", stringsAsFactors=F)
#opt <- data.frame(sample="m1035-108-dia", stringsAsFactors=F)
#opt <- data.frame(sample="m1060-108-rel", stringsAsFactors=F)
doplot <- FALSE
if (gtools::invalid(opt$sample)) stop("ERROR: sample not specified")

# determine sex
sex.file <- "/mnt/projects/p2ry8-crlf2/results/patient_sex.tsv"
sex.list <- read.delim(sex.file)
patient <- gsub("(C|D|R|R2|R3)$", "", opt$sample, perl=T)
sex <- sex.list$sex[sex.list$patient==patient]
if (gtools::invalid(sex)) stop(sprintf("ERROR: could not determine sex of patient %s (sample %s). Check if there is an entry in sex file %s", patient, opt$sample, sex.file))

# load counts and gene annotations
load("/mnt/projects/p2ry8-crlf2/results/exomeDepth/counts.RData")
load("/mnt/projects/p2ry8-crlf2/results/exomeCopy/genes.GRCh37v75.biomart.RData")
gr.genes <- GRanges(seqnames=genes$chromosome_name, ranges=IRanges(start=genes$start_position, end=genes$end_position), names=genes$hgnc_symbol)

# chromosome sizes
chrsize <- read.delim("/mnt/projects/generic/data/hg19/ucsc.hg19.chrom.sizes", check.names=FALSE, header=FALSE, col.names=c("chr", "size"), stringsAsFactors = FALSE)
chrsize$chr <- gsub("chr", "", chrsize$chr)

if (doplot) {
  pdf(paste0("/mnt/projects/p2ry8-crlf2/results/exomeDepth/", opt$sample, ".exomeDepth.pdf.part"), width=15, height=7)
}

cnvs <- NA
for (iteration in c("autosomes", "chr21", "XY")) {
  
  # convert to data frame, rename columns
  counts.df <- as.data.frame(counts)
  names(counts.df)[-6:-1] <- names(counts@values@unlistData)[-1]
  
  # filter exons depending on whether autosomes or sex chromosomes are processed
  if (iteration == "chr21") {
    counts.df <- counts.df[counts.df$space %in% c("21"),]
  } else if (iteration == "XY") {
    counts.df <- counts.df[counts.df$space %in% c("X", "Y"),]
  } else {
    counts.df <- counts.df[!counts.df$space %in% c("21", "X", "Y"),]
  }
  
  # select reference set
  test.counts <- counts.df[,opt$sample]
  ref.samples <- grep("C$", names(counts.df)[-6:-1], value = T)
  ref.samples <- ref.samples[ref.samples != opt$sample] # remove test sample from reference set if present
  if (iteration == "chr21") { # exclude trisomy 21 cases
    ref.samples <- ref.samples[!ref.samples %in% c("365C", "400C", "HW11537C", "887C", "360C", "506C", "1089C", "802C", "961C", "957C", "DL2C", "N7C", "DS10898C", "VS14645C", "SE15285C", "AL9890C", "GL11356C", "GI13C", "HV57C", "564C")]
  } else if (iteration == "XY") {  # only matched sex
    ref.samples <- ref.samples[gsub("(C|D|R|R2|R3)$", "", ref.samples, perl=T) %in% sex.list$patient[sex.list$sex==sex]]
  }
  ref.counts <- as.matrix(counts.df[,ref.samples])
  ref.choice <- select.reference.set(test.counts=test.counts, reference.counts=ref.counts, bin.length = counts.df$width, n.bins.reduced=10000)
  print(paste("Selected reference samples for sample", opt$sample, "iteration", iteration, ":", paste(ref.choice[[1]], collapse=",")))
  ref.matrix <- as.matrix(counts.df[,ref.choice$reference.choice, drop=FALSE])
  ref.selected <- apply(ref.matrix, MAR=1, FUN=sum)
  
  # call CNVs
  test.cnv <- new('ExomeDepth', test=test.counts, reference=ref.selected, formula='cbind(test, reference)~1')
  test.cnv <- CallCNVs(x=test.cnv, transition.probability = 10^-5, chromosome=counts.df$space, start=counts.df$start, end=counts.df$end, name=counts.df$names)
  #head(test.cnv@CNV.calls)
  
  if (nrow(test.cnv@CNV.calls) > 0) {
    # annotate common CNVs (Conrad et al., 2010)
    data(Conrad.hg19)
    test.cnv <- AnnotateExtra(x = test.cnv, reference.annotation = Conrad.hg19.common.CNVs, min.overlap = 0.5, column.name = 'Conrad.hg19')
  
    # annotate segments with overlapping gene names
    gr.cnvs <- makeGRangesFromDataFrame(test.cnv@CNV.calls, keep.extra.columns=TRUE, ignore.strand=TRUE, seqinfo=NULL, seqnames.field=c("chromosome"), start.field=c("start"), end.field=c("end"))
    gr.cnvs$genes <- as.character(NA)
    o <- findOverlaps(gr.cnvs, gr.genes)
    for(i in 1:nrow(test.cnv@CNV.calls)) { gnames <- values(gr.genes[o@subjectHits[o@queryHits==i]])$names; gr.cnvs$genes[i] <- paste(gnames[gnames!=""], collapse=",") }
    
    # add results of this iteration to combined results
    if (is.na(cnvs)) {
      cnvs <- as.data.frame(gr.cnvs)
    } else {
      cnvs <- rbind(cnvs, as.data.frame(gr.cnvs))
    }
  }

  # plot
  if (doplot) {
    for (chr in unique(test.cnv@CNV.calls$chromosome)) {
      if (chr == "MT") { next }
      plot(test.cnv, sequence=chr, xlim=c(1, chrsize$size[chrsize$chr==chr]), count.threshold=20, main=paste(opt$sample, "chromosome", chr), cex.lab=1, with.gene=FALSE)
    }
  }
}

if (doplot) {
  dev.off()
}

# write results
cnvs$sample.name <- opt$sample
cnvs <- cnvs[order(cnvs$BF, decreasing=T),]
write.table(cnvs, file=paste0("/mnt/projects/p2ry8-crlf2/results/exomeDepth/", opt$sample, ".cnvs.tsv.part"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

