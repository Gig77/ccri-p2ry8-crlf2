options(warn=1)
library(GenomicRanges)
library(reshape)

files <- list.files(path="/mnt/projects/p2ry8-crlf2/results/exomeDepth", pattern="*.cnvs.tsv$", full.names=T)

d <- read.delim(files[1], stringsAsFactor=F)
for (i in 2:length(files)) {
	d <- rbind(d, read.delim(files[i], stringsAsFactor=F))
}

d <- d[,!names(d) %in% c("strand", "width")]
gr <- makeGRangesFromDataFrame(d, keep.extra.columns=TRUE, ignore.strand=TRUE, seqinfo=NULL, seqnames.field=c("seqnames"), start.field=c("start"), end.field=c("end"))

# find overlaps
o <- findOverlaps(gr, gr)
o <- o[o@queryHits != o@subjectHits] # remove self-hits
#o <- o[gr$copy.count[o@queryHits] != 2 & gr$copy.count[o@subjectHits] != 2] # remove segments with normal copy number two
o <- o[gr$type[o@queryHits] == gr$type[o@subjectHits]] 
o <- o[width(gr)[o@queryHits] >= 20 & width(gr)[o@subjectHits] >= 20] # remove overlaps with tiny single-exon segments

# ignore overlaps with large segments in normal samples due to constitutional trisomy 21 
#const.tri21 <- which(seqnames(gr)=="21" & grepl("C$", gr$sample.name, perl=T) & gr$type == "duplication" & gr@ranges@width >= 1000000)
#o <- o[!(o@queryHits %in% const.tri21) & !(o@subjectHits %in% const.tri21)]

# remove overlaps with samples of questionable quality (= correlation with reference sets < 0.95)
crappy <- c("m1059-92-dia", "B36D", "LU3D", "B36R", "m1037-839-dia", "92R", "92D", "GI8R", "m1037-839-dia", "LU3R", "m247-833-dia", "m248-841-dia", "m1060-108-rel", 
            "MA5R", "S23R3", "360D", "m1069-737-rel", "1089D", "m1035-108-dia", "839D", "S23R", "m1069-737-rel", "m1977-G-dia", "108C", "839R", "m1964-545-rel", "m252-379-dia")
#crappy <- c("MA5R", "S23R", "242C", "HV57R", "S23R3", "KE17247R")
o <- o[!(gr$sample.name[o@subjectHits] %in% crappy),]
#o <- o[!(seqnames(gr)[o@subjectHits] %in% c("X","Y") & gr$sample.name[o@subjectHits] %in% c("GI8C")),] # GI8C gives weird results for sex chromosomes...

# determine overlap in percent of shared exons
ex <- read.delim("/mnt/projects/generic/data/illumina/nexterarapidcapture_exome_targetedregions.nochr.bed", header=F)
er <- GRanges(seqnames=ex$V1, ranges=IRanges(start=ex$V2, end=ex$V3))
oex <- findOverlaps(gr, er)
segex <- cast(as.data.frame(oex), formula=queryHits~., value="subjectHits", fun.aggregate=function(x) { paste(x, collapse=",") })
names(segex) <- c("subjectHits", "subjectExons")
om <- merge(as.data.frame(o), segex)
names(segex) <- c("queryHits", "queryExons")
om <- merge(om, segex)
om <- om[order(om$queryHits, om$subjectHits),]
o@metadata[["pct_shared_exons"]] <- apply(om, 1, function(x) { qex <- unlist(strsplit(x["queryExons"], ",", fixed=T)) ; sex <- unlist(strsplit(x["subjectExons"], ",", fixed=T)) ; length(intersect(qex, sex)) / length(union(qex, sex)) })
o <- o[o@metadata[["pct_shared_exons"]] >= 0.3]

# to be considered same event, we require overlap of at least 30% (intersection divided by union)
#qr <- gr[o@queryHits]@ranges
#sr <- gr[o@subjectHits]@ranges
#o@metadata[["width"]] <- (pmin(end(qr), end(sr))-pmax(start(qr), start(sr))+1) / (pmax(end(qr), end(sr))-pmin(start(qr), start(sr))+1)
#o@metadata[["pct_shared_exons"]] <- (pmin(end(qr), end(sr))-pmax(start(qr), start(sr))+1) / (pmax(end(qr), end(sr))-pmin(start(qr), start(sr))+1)
#o <- o[o@metadata[["width"]] >= 0.3]

gr$overlap.samples <- as.character(NA)
gr$overlap.count <- as.integer(NA)
gr$overlap.count.tumor <- as.integer(NA)
library(plyr)
o.aggregated <- ddply(as.data.frame(o), .(queryHits), function(x) { 
			c(overlap.samples=paste(unique(gr$sample.name[x$subjectHits]), collapse=","), 
			  overlap.count=length(unique(gr$sample.name[x$subjectHits])),
			  overlap.count.tumor=sum(!grepl("C$", unique(gr$sample.name[x$subjectHits])))) 
  })
gr$overlap.samples[o.aggregated$queryHits] <- o.aggregated$overlap.samples
gr$overlap.count[o.aggregated$queryHits] <- o.aggregated$overlap.count
gr$overlap.count.tumor[o.aggregated$queryHits] <- o.aggregated$overlap.count.tumor

# write table
write.table(as.data.frame(gr), file="/mnt/projects/p2ry8-crlf2/results/exomeDepth/allsamples.cnvs.exomeDepth.tsv.part", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, na="")

