library(Gviz)
library(GenomicRanges)

#region <- "IKZF1"; chr <- "chr7"; start <- 50000000; end <- 50800000
#region <- "CDKN2AB"; chr <- "chr9"; start <- 21600000; end <- 22400000
#region <- "PAR1"; chr <- "chrX"; start <- 1000000; end <- 2000000
#region <- "PAX5"; chr <- "chr9"; start <- 36495626; end <- 37377382
region <- "GABRB3"; chr <- "chr15"; start <- 25870574; end <- 27936343
#region <- "ACSM2A"; chr <- "chr16"; start <- 20318327; end <- 20643523
#region <- "BTG1"; chr <- "chr12"; start <- 91278699; end <- 92889529

s <- read.delim("~/p2ry8-crlf2/results/exomeCopy/allpatients.filtered-segments.exomeCopy.tsv")
gr <- GRanges(seqnames=paste0("chr", s$seqnames), ranges=IRanges(start=s$start, end=s$end), sample.name=s$sample.name, group=s$sample.name, copy.count=s$copy.count, log.odds=s$log.odds, nranges=s$nranges, targeted.bp=s$targeted.bp, genes=s$genes)
gr <- gr[seqnames(gr)==chr & end(gr) >= start & start(gr) <= end,]
start(gr) <- pmax(start(gr), start)
end(gr) <- pmin(end(gr), end)

# get gene models in region of interest
options(ucscChromosomeNames=FALSE)
biomTrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = chr, start = start, end = end, name = "Genes", showId=TRUE)

# keep only longest transcript per gene
trlen <- aggregate(width ~ gene + transcript, data=as.data.frame(biomTrack@range), FUN=sum)
biomTrack@range <- biomTrack@range[biomTrack@range$transcript %in% trlen[ave(trlen$width, trlen$gene, FUN=max)==trlen$width,"transcript"]]
seqlevels(ranges(biomTrack)) <- paste0("chr", seqlevels(ranges(biomTrack)))

# plot
itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
gtrack <- GenomeAxisTrack()

cntracks <- list()
for (sample in unique(gr$sample.name)) {
	print(sample)
	grs <- gr[gr$sample.name==sample,]
	cntrack <- AnnotationTrack(grs, name=sample, stackHeight=1, rot.title=0, cex.title=0.5, col.line="white")
	feature(cntrack) <- paste0("CN", grs$copy.count) 
	cntracks <- c(cntracks, cntrack)
}

pdf(paste0("~/p2ry8-crlf2/results/exomeCopy/", region, ".pdf"))
plotTracks(c(list(itrack, gtrack, biomTrack), cntracks), from=start, to=end, title.width=1.8, CN0="darkred", "CN1"="red", "CN2"="white", "CN3"="lightblue", "CN4"="darkblue")
dev.off()

#cntrack <- AnnotationTrack(gr, name="CN", showId=T, stacking="squish", col="black", cex=0.7, stackHeight=1, fontcolor="black", just.group="above", rot.title=0)
#cntrack@range@elementMetadata@listData$id <- gr$sample.name

#pdf("~/p2ry8-crlf2/results/exomeCopy/IKZF1.pdf")
#plotTracks(list(itrack, gtrack, biomTrack, cntrack), from=start, to=end)
#dev.off()
