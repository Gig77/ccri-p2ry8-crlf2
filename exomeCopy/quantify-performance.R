options(warn=1)
library(GenomicRanges)

exons.df <- read.delim("~/generic/data/illumina/nexterarapidcapture_exome_targetedregions.nochr.bed", header = FALSE)
exons <- GRanges(seqname = exons.df[, 1], IRanges(start = exons.df[,2] + 1, end = exons.df[, 3]))

snp.df <- read.delim("~/p2ry8-crlf2/results/exomeCopy/P2RY8-CRLF2_SNParray.txt")
snp.df <- snp.df[snp.df$Type != "LOH",]
snp <- GRanges(seqname = snp.df[,"Chr"], IRanges(start = snp.df[,"Min"], end = snp.df[,"Max"]), cn = snp.df[,"CN.State"], sample = snp.df[,"Sample"])

ec.df <- read.delim("~/p2ry8-crlf2/results/exomeCopy/allpatients.filtered-segments.exomeCopy.tsv")
ec <- GRanges(seqname = ec.df[,"seqnames"], IRanges(start = ec.df[,"start"], end = ec.df[,"end"]), cn = ec.df[,"copy.count"], sample = ec.df[,"sample.name"])

s <- "737R"
name.snp <- paste0(s,"a") ; name.ec <- paste0(s,"e")

snp.sample <- snp[snp$sample==s,]
ec.sample <- ec[ec$sample==s,]
o.snp <- findOverlaps(exons, snp.sample)
o.ec <- findOverlaps(exons, ec.sample)

mcols(exons)[,name.snp] <- 2
mcols(exons[o.snp@queryHits])[,name.snp] <- snp.sample[o.snp@subjectHits]$cn
mcols(exons)[,name.ec] <- 2
mcols(exons[o.ec@queryHits])[,name.ec] <- ec.sample[o.ec@subjectHits]$cn

sens <- sum(exons[, name.snp] != 2 & mcols(exons)[, name.ec] != 2) / sum(mcols(exons)[, name.snp] != 2)
print(sprintf("Sensitivity: %.2f%%", sens*100))
spec <- sum(mcols(exons)[, name.snp] == 2 & mcols(exons)[, name.ec] == 2) / sum(mcols(exons)[, name.snp] == 2)
print(sprintf("Specificity: %.2f%%", spec*100))
prec <- sum(mcols(exons)[, name.snp] != 2 & mcols(exons)[, name.ec] != 2) / sum(mcols(exons)[, name.ec] != 2)
print(sprintf("Precision or positive predictive value (PPV): %.2f%%", prec*100))
