options(warn=1)
library(optparse)
library(gtools)

option_list <- list(
			make_option("--segments-file", type="character", help="file with exomeCopy segments"),
			make_option("--output-file", type="character", help="output file name"),
			make_option("--tumor", type="character", help="sample name tumor"),
			make_option("--normal", type="character", help="sample name normal")
		)

opt <- parse_args(OptionParser(option_list=option_list))
#opt <- data.frame('segments-file'="/mnt/projects/p2ry8-crlf2/results/exomeCopy/allpatients.compiled-segments.exomeCopy.tsv", 'output-file'="/mnt/projects/p2ry8-crlf2/results/exomeCopy/circos/S737R2.somatic.circos.tsv", tumor="S737R2", normal="S737C", stringsAsFactors=F, check.names=F)

if (invalid(opt$'segments-file')) stop("segments file not specified")
if (invalid(opt$'output-file')) stop("output file not specified")
if (invalid(opt$tumor)) stop("tumor sample name not specified")
if (invalid(opt$normal)) stop("normal sample name not specified")

ca <- read.delim(opt$'segments-file', stringsAsFactor=F)
ca <- ca[ca$sample==opt$tumor,]

if (nrow(ca) > 0) {
	colors <- list("0"="dred", "1"="red", "2"="black", "3"="blue", "4"="dblue", "5"="vdblue", "6"="vvdblue")
	out <- data.frame(paste0("hs", ca$seqnames), ca$start, ca$end, paste0("fill_color=", sapply(colors[as.character(ca$copy.count)], "[[", 1)))
	write.table(out, file=opt$'output-file', quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)	
} else {
	system(paste("touch", opt$'output-file')) # create empty file
}
