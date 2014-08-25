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
#opt <- data.frame('segments-file'="~/p2ry8-crlf2/results/exomeCopy/allpatients.compiled-segments.exomeCopy.tsv", 'output-file'="~/p2ry8-crlf2/results/exomeCopy/circos/S737R2.somatic.circos.tsv", tumor="S737R2", normal="S737C", stringsAsFactors=F, check.names=F)

if (invalid(opt$'segments-file')) stop("segments file not specified")
if (invalid(opt$'output-file')) stop("output file not specified")
if (invalid(opt$tumor)) stop("tumor sample name not specified")
if (invalid(opt$normal)) stop("normal sample name not specified")

ca <- read.delim(opt$'segments-file', stringsAsFactor=F)
cf <- ca[ca$sample==opt$tumor & ca$copy.count != 2 & (ca$log.odds == 0 | ca$log.odds > 20) & !grepl("C", ca$overlap.samples),]

if (nrow(cf) > 0) {
	colors <- list("0"="dred", "1"="red", "2"="black", "3"="blue", "4"="dblue", "5"="vdblue", "6"="vvdblue")
	out <- data.frame(paste0("hs", cf$seqnames), cf$start, cf$end, paste0("fill_color=", sapply(colors[as.character(cf$copy.count)], "[[", 1)))
	write.table(out, file=opt$'output-file', quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)	
} else {
	system(paste("touch", opt$'output-file')) # create empty file
}
