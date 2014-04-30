options(warn=1)
library(optparse)
library(gtools)
library(DNAcopy)

set.seed(25)
min.cov.rem <- 10

option_list <- list(
		make_option("--patient", type="character", help="patient ID"),
		make_option("--diagnosis", type="character", help="diagnosis coverage data file"),
		make_option("--relapse", type="character", help="relapse coverage data file"),
		make_option("--remission", type="character", help="remission coverage data file"),
		make_option("--output", type="character", help="PDF output file"),
		make_option("--region-name", type="character", help="region name to display"),
		make_option("--display-genes", type="character", help="gene names to display"),
		make_option("--display-chrom", type="character", help="chromosome to display"),
		make_option("--display-start", type="character", help="start coordinate of region to display"),
		make_option("--display-end", type="character", help="end coordinate of region to display")
		)
opt <- parse_args(OptionParser(option_list=option_list))
#opt <- data.frame(patient="92", diagnosis="coverage/92D.exon-coverage.tsv", relapse="coverage/92R.exon-coverage.tsv", remission="coverage/92C.exon-coverage.tsv", output="coverage/92.PAR.pdf", stringsAsFactors=F)

if (invalid(opt$patient)) stop("patient not specified")
if (invalid(opt$diagnosis) && invalid(opt$relapse)) stop("diagnosis sample, relapse sample, or both need to be specified")
if (invalid(opt$remission)) stop("remission sample not specified")
if (invalid(opt$output)) stop("output file not specified")
if (invalid(opt$'region-name')) stop("region name not specified")
if (invalid(opt$'display-chrom')) stop("chromosome not specified")
if (invalid(opt$'display-start')) stop("start coordinate not specified")
if (invalid(opt$'display-end')) stop("end coordinate not specified")

region.chr <- opt$'display-chrom' # e.g. "X"
region.start <- as.numeric(opt$'display-start') # e.g. 60001
region.end <- as.numeric(opt$'display-end') # e.g. 2699520

display.genes <- c()
if (!invalid(opt$'display-genes')) display.genes <- strsplit(opt$'display-genes', ",")[[1]]

rem <- read.delim(opt$remission, check.names=F, stringsAsFactor=F, header=F)
names(rem) <- c("chr", "start", "end", "name", "length", "strand", "total.rem", "avg.rem")
avg.rem <- mean(rem$avg.rem)
rem <- rem[rem$chr==region.chr & rem$start >= region.start & rem$end <= region.end,]

dia <- read.delim(opt$diagnosis, check.names=F, stringsAsFactor=F, header=F)
names(dia) <- c("chr", "start", "end", "name", "length", "strand", "total.dia", "avg.dia")
avg.dia <- mean(dia$avg.dia)
dia <- dia[dia$chr==region.chr & dia$start >= region.start & dia$end <= region.end,]
m <- merge(rem, dia, by=c("chr", "start", "end", "name", "length", "strand"))
m$ratio.dia <- log2(m$avg.dia/m$avg.rem)-log2(avg.dia/avg.rem)

# segment data
CNA.object.dia <-CNA(genomdat = m$ratio.dia, chrom = m$chr, maploc = m$start, data.type = 'logratio')
CNA.smoothed.dia <- smooth.CNA(CNA.object.dia)
#segs <- segment(CNA.smoothed, alpha = 0.01, verbose=0, min.width=2, undo.splits="sdundo", undo.SD=3)
segs.dia <- segment(CNA.smoothed.dia, alpha = 0.01, verbose=0, min.width=2)

if (!invalid(opt$relapse)) {
	rel <- read.delim(opt$relapse, check.names=F, stringsAsFactor=F, header=F)
	names(rel) <- c("chr", "start", "end", "name", "length", "strand", "total.rel", "avg.rel")
	avg.rel <- mean(rel$avg.rel)
	rel <- rel[rel$chr==region.chr & rel$start >= region.start & rel$end <= region.end,]
	m <- merge(m, rel, by=c("chr", "start", "end", "name", "length", "strand"))
	m$ratio.rel <- log2(m$avg.rel/m$avg.rem)-log2(avg.rel/avg.rem)

	# segment data
	CNA.object.rel <-CNA(genomdat = m$ratio.rel, chrom = m$chr, maploc = m$start, data.type = 'logratio')
	CNA.smoothed.rel <- smooth.CNA(CNA.object.rel)
	#segs <- segment(CNA.smoothed, alpha = 0.01, verbose=0, min.width=2, undo.splits="sdundo", undo.SD=3)
	segs.rel <- segment(CNA.smoothed.rel, alpha = 0.01, verbose=0, min.width=2)
}

# remove exons below minimum coverage in remission
#m <- m[m$total.rem/m$length>=min.cov.rem,]

# get gene coordinates
dia$gene <- sapply(strsplit(as.character(dia$name), ":"), "[", 3)
genes <- data.frame(gene=aggregate(start~gene, dia, min)$gene, start=aggregate(start~gene, dia, min)$start, end=aggregate(end~gene, dia, max)$end)
if (length(display.genes) > 0) genes <- genes[genes$gene %in% display.genes,]
genes <- genes[order(genes$start),]

# normalize coverage to total coverage
mn <- m
#mn <- mn[mn$avg.dia >= 20 | mn$avg.rel >= 20 | mn$avg.rem >= 20,]
#mn$avg.dia <- mn$avg.dia / median(mn[mn$chr=="chr19", "avg.dia"]) * median(mn[mn$chr=="chr19", "avg.rem"])
#mn$avg.rel <- mn$avg.rel / median(mn[mn$chr=="chr19", "avg.rel"]) * median(mn[mn$chr=="chr19", "avg.rem"])

# plot pseudoautosomal region 1 (PAR1)
pdf(opt$output, height=8)
par(mfrow=c(2,1), mar=c(4.5,4,3,2))
plot(mn$start, mn$ratio.dia, ylim=c(-2, 2), xlim=c(region.start,region.end), col=rgb(0,0,0,1), main=paste(opt$patient, opt$'region-name', "dia"), ylab="log2 cov ratio", xlab=paste0(region.chr, ":", region.start, "-", region.end), cex=0.5)
for(i in 1:length(segs.dia$output$loc.start)) {
	lines(c(segs.dia$output$loc.start[i], segs.dia$output$loc.end[i]),c(segs.dia$output$seg.mean[i],segs.dia$output$seg.mean[i]),type="l", col="orange", lty=1, lwd=3)
}
text(x=(genes$start+genes$end)/2, y=rep(c(-1.8, -1.95), length(genes$start)/2), labels=genes$gene, cex=0.6)

if (!invalid(opt$relapse)) {
	plot(mn$start, mn$ratio.rel, ylim=c(-2, 2), xlim=c(region.start,region.end), col=rgb(0,0,0,1), main=paste(opt$patient, opt$'region-name', "rel"), ylab="log2 cov ratio", xlab=paste0(region.chr, ":", region.start, "-", region.end), cex=0.5)
	for(i in 1:length(segs.rel$output$loc.start)) {
		lines(c(segs.rel$output$loc.start[i], segs.rel$output$loc.end[i]),c(segs.rel$output$seg.mean[i],segs.rel$output$seg.mean[i]),type="l", col="orange", lty=1, lwd=3)
	}
	text(x=(genes$start+genes$end)/2, y=rep(c(-1.8, -1.95), length(genes$start)/2), labels=genes$gene, cex=0.6)
}

dev.off()

