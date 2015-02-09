library(lattice)
library(optparse)
library(DNAcopy)

option_list <- list(make_option("--input-file", type="character", help="name of input file"),
					make_option("--output-file", type="character", help="name of output file"),
					make_option("--input-file-normal", type="character", help="name of input file (matched normal)"))
opt <- parse_args(OptionParser(option_list=option_list))
			
#opt <- data.frame('input-file' = "~/p2ry8-crlf2/results/varscan/715R.varscan.dbsnp.vcf", 'input-file-normal' = "~/p2ry8-crlf2/results/varscan/715C.varscan.dbsnp.vcf", 'output-file' = "~/p2ry8-crlf2/results/loh/715R.loh.pdf", stringsAsFactors=F, check.names=F)	
			
v <- read.table(opt$'input-file', header=F, sep="\t", stringsAsFactors=F)
v <- v[v[,1] %in% c(1:22, "X", "Y") & grepl("rs", v[,3]),]

dp <- unlist(lapply(strsplit(v[,10], ":"), "[", 4))
pval <- as.numeric(unlist(lapply(strsplit(v[,10], ":"), "[", 8)))
af <- as.numeric(gsub("\\%", "", unlist(lapply(strsplit(v[,10], ":"), "[", 7)))) / 100

d <- data.frame(chr=factor(v[,1], levels=c(1:22, "X", "Y")), pos=as.numeric(v[,2]) / 1000000, af, dbsnp=v[,3])
d <- d[pval<1e-6,]

# remove homoyzgous SNPs if matched normal is provided
if (!is.null(opt$'input-file-normal')) {
	vn <- read.table(opt$'input-file-normal', header=F, sep="\t", stringsAsFactors=F)
	vn <- vn[vn[,1] %in% c(1:22, "X", "Y") & grepl("rs", vn[,3]),]
	afn <- as.numeric(gsub("\\%", "", unlist(lapply(strsplit(vn[,10], ":"), "[", 7)))) / 100
	pvaln <- as.numeric(unlist(lapply(strsplit(vn[,10], ":"), "[", 8)))
	homsnp <- vn[afn>=0.8,3]
	keep <- (d$dbsnp %in% vn$V3) & (!d$dbsnp %in% homsnp)
	print(sprintf("%d homozygous SNPs removed", sum(!keep)))
	d <- d[keep,]
}

print(sprintf("Number of SNPs: %d", nrow(d)))

# segmentation of "mirrored" BAF
d$maf <- ifelse(is.null(opt$'input-file-normal') & d$af > 0.9, NA, abs(d$af - 0.5) + 0.5)
CNA.object <-CNA(genomdat = d[,"maf"], chrom = d[,"chr"], maploc = d[,"pos"], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, alpha = 0.001, verbose=0, min.width=2)


pdf(opt$'output-file', width=30, paper="A4r")
#xyplot(af~pos | chr, d, cex=0.3, as.table=T, col="black", layout=c(2, 12), ylim=c(0, 1))
par(mfrow=c(5,5), mar=c(0.5,2,1.5,1), oma=c(0, 0, 2, 0))
for (c in levels(d$chr))
{
	if (length(d[d$chr==c,"pos"]) > 0) {
		plot(d[d$chr==c,"pos"], d[d$chr==c,"af"], ylim=c(0,1), main=c, cex=0.3, xaxt='n', yaxt='n', col="darkblue", lwd=0.5)
	}
	else {
		plot(0, 0, ylim=c(0,1), main=c, cex=0.3, xaxt='n', yaxt='n', col="darkblue", lwd=0.5)
	}
	for(i in which(segs$output$chrom==c & segs$output$seg.mean >= 0.6)) {
		lines(c(segs$output$loc.start[i], segs$output$loc.end[i]),c(segs$output$seg.mean[i],segs$output$seg.mean[i]),type="l", col="orange", lty=1, lwd=2)
		lines(c(segs$output$loc.start[i], segs$output$loc.end[i]),c(1-segs$output$seg.mean[i],1-segs$output$seg.mean[i]),type="l", col="orange", lty=1, lwd=2)
	}
	axis(2, at=c(0, 0.25, 1/3, 0.5, 2/3, 0.75, 1), labels=c("0", "", "", "0.5", "", "", "1"))
	abline(h=0.5)
	abline(h=c(1/3, 0.25, 2/3, 0.75), lty=2)
}
mtext(opt$'input-file', outer=TRUE, cex=1.5)
dev.off()
