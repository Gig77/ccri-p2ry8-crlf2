options(warn=1)
library(GenomicRanges)

exons.df <- read.delim("/mnt/projects/generic/data/illumina/nexterarapidcapture_exome_targetedregions.nochr.bed", header = FALSE)
exons <- GRanges(seqname = exons.df[, 1], IRanges(start = exons.df[,2] + 1, end = exons.df[, 3]))

snp.df <- read.delim("/mnt/projects/p2ry8-crlf2/results/exomeCopy/P2RY8-CRLF2_SNParray.txt")
snp.df <- snp.df[snp.df$Type != "LOH",]
snp <- GRanges(seqname = snp.df[,"Chr"], IRanges(start = snp.df[,"Min"], end = snp.df[,"Max"]), cn = snp.df[,"CN.State"], sample = snp.df[,"Sample"])

ec.df <- read.delim("/mnt/projects/p2ry8-crlf2/results/exomeCopy/allpatients.filtered-segments.exomeCopy.tsv")
ec <- GRanges(seqname = ec.df[,"seqnames"], IRanges(start = ec.df[,"start"], end = ec.df[,"end"]), cn = ec.df[,"copy.count"], sample = ec.df[,"sample.name"])

perf <- data.frame(sample=character(0), tp=integer(0), fn=integer(0), fp=integer(0), tn=integer(0), sens=numeric(0), spec=numeric(0), prec=numeric(0), stringsAsFactors=F)

for(s in c("242D", "379D", "737D", "737R", "839D", "839R", "92D", "92R", "957D", "961D")) {
	name.snp <- paste0(s,"a")
	name.ec <- paste0(s,"e")
	
	snp.sample <- snp[snp$sample==s,]
	ec.sample <- ec[ec$sample==s,]
	o.snp <- findOverlaps(exons, snp.sample)
	o.ec <- findOverlaps(exons, ec.sample)
	
	mcols(exons)[,name.snp] <- 2
	mcols(exons[o.snp@queryHits])[,name.snp] <- snp.sample[o.snp@subjectHits]$cn
	mcols(exons)[,name.ec] <- 2
	mcols(exons[o.ec@queryHits])[,name.ec] <- ec.sample[o.ec@subjectHits]$cn
	
	tp <- sum(mcols(exons)[, name.snp] != 2 & mcols(exons)[, name.ec] != 2)
	fn <- sum(mcols(exons)[, name.snp] != 2 & mcols(exons)[, name.ec] == 2)
	fp <- sum(mcols(exons)[, name.snp] == 2 & mcols(exons)[, name.ec] != 2)
	tn <- sum(mcols(exons)[, name.snp] == 2 & mcols(exons)[, name.ec] == 2)
	
	sens <- tp / (tp + fn)
	spec <- tn / (tn + fp)
	prec <- tp / (tp + fp)
	
	perf[nrow(perf)+1,] <- c(s, tp, fn, fp, tn, sens, spec, prec)
}

write.table(perf, "/mnt/projects/p2ry8-crlf2/results/exomeCopy/performance.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
