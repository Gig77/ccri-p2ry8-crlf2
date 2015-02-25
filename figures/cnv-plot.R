min.width <- 3000000
# chromosome sizes
sizes <- read.delim("~/generic/data/hg19/ucsc.hg19.chrom.sizes", header=F, colClasses=c("factor", "numeric"))
colnames(sizes) <- c("chr", "size")
sizes$chr <- gsub("chr", "", sizes$chr)
rownames(sizes) <- sizes$chr
sizes <- sizes[c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),]
sizes$offset <- cumsum(sizes$size) - sizes$size

# sex
s <- read.delim("~/p2ry8-crlf2/results/patient_sex.tsv")

# cohorts
samples.matched <- c("DL2D", "DL2R", "108D", "108R", "108R2", "GI8D", "GI8R", "HV80D", "HV80R","N7D", "N7R", "DS10898D", "DS10898R", "VS14645D", "VS14645R", "BB16D", "BB16R", "SE15285D", "SE15285R", "715D", "715R", "715R3", "839D", "839R", "AL9890D", "AL9890R", "GL11356D", "GL11356R", "GI13D", "GI13R", "B36D", "B36R", "92D", "92R", "HV57D", "HV57R", "737D", "737R", "737R3", "KE17247D", "KE17247R", "S23D", "S23R", "1060D", "BJ17183D")
cases.matched <- unique(gsub("(D|R\\d?)$", "", samples.matched))
samples.diaonly <- c("242D", "360D", "365D", "379D", "400D", "506D", "769D", "802D", "833D", "887D", "841D", "903D", "948D", "961D", "1060D", "1066D", "1089D", "HW11537D", "KT14158D", "TL14516D")
samples.mlpa <- c("1060D", "1066D", "1089D", "108D", "108R", "242D", "360D", "365D", "460D", "460R", "545D", "545R", "737R", "769D", "802D", "833D", "839D", "887D", "903D", "948D", "957D", "961D", "B36R", "BB16D", "BB16R", "DL2D", "DL2R", "DS10898D", "DS10898R", "GI13R", "GI8D", "HV57D", "HV57R", "HV80D", "HV80R", "KE17247R", "KT14158D", "SE15285D", "SE15285R", "SN18D", "SN18R", "VS14645D", "VS14645R", "564D", "715D", "715R", "GI8R")
cases.excluded <- c("LU3", "SN18", "564", "460", "545", "MA5", "957")
cases.ds <- c("DL2", "N7", "DS10898", "VS14645", "SE15285", "AL9890", "GL11356", "GI13", "HV57", "564", "365", "400", "HW11537", "887", "360", "506", "1089", "802", "961")
samples.lowqual <- c("GI8R", "HV57R", "KE17247R", "BJ17183R", "1060R")

# marker genes
markers <- data.frame(chr="7", pos=50344378, name="IKZF1", padj=0.5, stringsAsFactors=F)
markers <- rbind(markers, c("9", 21967751, "CDKN2A/B", 0.1))
markers <- rbind(markers, c("6", 26250835, "HIST1", 0.5))
markers <- rbind(markers, c("9", 36966544, "PAX5", 0.9))
markers <- rbind(markers, c("X", 1581466, "PAR1", 0.5))
markers <- rbind(markers, c("12", 92534054, "BTG1", 0.5))
markers <- rbind(markers, c("4", 87967845, "AFF1", 0.5))
#markers <- rbind(markers, c("3", 47057898, "SETD2", 0.5))
markers <- rbind(markers, c("5", 158139160, "EBF1", 0.5))
markers$offset <- as.numeric(markers$pos) + sizes[match(markers$chr, sizes$chr), "offset"]

# read cnvs
d <- read.delim("~/p2ry8-crlf2/results/exomeCopy/allpatients.filtered-segments.exomeCopy.tsv", stringsAsFactors=F)
d.hdall <- read.delim("~/hdall/results/exomeCopy/allpatients.filtered-segments.exomeCopy.tsv", stringsAsFactors=F)
d.hdall$seqnames <- gsub("chr", "", d.hdall$seqnames)
d.hdall$sample.name <- gsub("_dia", "D", d.hdall$sample.name)
d.hdall$sample.name <- gsub("_rel", "R", d.hdall$sample.name)
d <- rbind(d, d.hdall)
d$source <- "WES"

# add SNParray results
d.array <- read.delim("~/p2ry8-crlf2/results/cnvs.snp_arrays.txt", stringsAsFactor=F)
d <- rbind(d, data.frame(source="array", seqnames=d.array$seqnames, start=d.array$start, end=d.array$end, width=d.array$size, strand="*", sample.name=d.array$sample, copy.count=d.array$copy.count, log.odds=NA, nranges=d.array$marker, targeted.bp=NA, genes=d.array$genes, overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA))

# exclude IGH/TCR regions
d <- d[!(d$seqnames == "7" & d$start >= 38000000 & d$end <= 38500000),] # TCR
d <- d[!(d$seqnames == "7" & d$start >= 141000000 & d$end <= 143000000),] # TCR
d <- d[!(d$seqnames == "14" & d$start >= 22000000 & d$end <= 24000000),] # TCR
#d <- d[!(d$seqnames == "22" & d$start >= 22000000 & d$end <= 24000000),] # IGLV
d <- d[!(d$seqnames == "2" & d$start >= 88000000 & d$end <= 91000000),] # IGK
d <- d[!(d$seqnames == "2" & d$start >= 97600000 & d$end <= 98300000),] # IGK
d <- d[!(d$seqnames == "14" & d$start >= 102000000 & d$end <= 108000000),] # IGH

# exclude germline trisomy
d <- d[!(d$sample.name %in% c(paste0(cases.ds, "D"), paste0(cases.ds, "R"))) | d$seqnames != "21" | d$copy.count != 3,]

# add MLPA CNVs missed by exomeCopy
d <- rbind(d,data.frame(source="MLPA", seqnames="7", start=50435703, end=50459561, width=23859, strand="*", sample.name="108R", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="IKZF1", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="7", start=50435703, end=50459561, width=23859, strand="*", sample.name="DL2R", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="IKZF1", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="7", start=50358657, end=50459561, width=100905, strand="*", sample.name="DS10898D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="IKZF1", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="7", start=50358657, end=50459561, width=100905, strand="*", sample.name="GI8R", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="IKZF1", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="9", start=21968226, end=22029516, width=61291, strand="*", sample.name="360D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="CDKN2A", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="9", start=21968226, end=22029516, width=61291, strand="*", sample.name="957D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="CDKN2A", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="9", start=21968226, end=22029516, width=61291, strand="*", sample.name="SN18D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="CDKN2A", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="9", start=21968226, end=22029516, width=61291, strand="*", sample.name="SN18R", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="CDKN2A", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="9", start=21968226, end=22029516, width=61291, strand="*", sample.name="GI8R", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="CDKN2A", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="9", start=36923350, end=37026572, width=103223, strand="*", sample.name="564D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="PAX5", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="9", start=36923350, end=37026572, width=103223, strand="*", sample.name="564R", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="PAX5", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="GI8D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="GI8R", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="961D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="715D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="715R", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="MLPA", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="GL11356D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 

# add cryptic P/C with copy number two b/c of X gain
d <- rbind(d,data.frame(source="WES", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="715R3", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="WES", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="108D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="WES", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="108R", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="WES", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="108R2", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="WES", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="N7R", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="WES", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="506D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="WES", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="903D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
d <- rbind(d,data.frame(source="WES", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="1066D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 

# add WES false-negative events based on visual inspection
d <- rbind(d,data.frame(source="WES", seqnames="X", start=1387692, end=1713119, width=325428, strand="*", sample.name="802D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="P2RY8,CRLF2", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 

# add cumulative genome start/end coordinate
d$start.cum <- d$start + sizes[match(d$seqnames, sizes$chr), "offset"]
d$end.cum <- d$end + sizes[match(d$seqnames, sizes$chr), "offset"]

# assign colors
col.0 <- "black"
col.1 <- rgb(57/255,76/255,155/255)
col.3 <- rgb(237/255,133/255,151/255)
col.4 <- rgb(232/255,44/255,40/255)
col.5 <- "orange"
d$color <- NA
d$color[d$copy.count == 0] <- col.0
d$color[d$copy.count == 1] <- col.1
d$color[d$copy.count == 3] <- col.3
d$color[d$copy.count == 4] <- col.4
d$color[d$copy.count >= 5] <- col.5
d$color[d$seqnames %in% c("X", "Y") & d$copy.count == 2 & s[match(gsub("(.*)(D|R\\d?)$", "\\1", d$sample.name), s$patient), "sex"] %in% c("m")] <- col.3 # add. sex chr. in males
d$color[d$seqnames %in% c("X", "Y") & d$copy.count > 2] <- col.3 # smooth exomeCopy calls oscillating between 3 and 4 copies 

plot.cnvs <- function(samples) {
	plot(0,0,xlim=c(0,sum(sizes$size)),ylim=c(1,numsamp+1),type='n',xaxt='n',yaxt='n',xlab="", ylab="", xaxs="i", yaxs="i")

	# gray background for samples with missing data, axis to the right
	for (i in c(1:numsamp)) {	
		sample <- samples[i]
		
		# annotate CNV source on right y axis
		source <- NULL
		if (sample %in% c(unique(d.array$sample))) {
			source <- "A"
		} else if (!sample %in% c(samples.lowqual)) {
			source <- "E"
		} else {
			rect(0, (numsamp-i+1)+0.98, sum(sizes$size), (numsamp-i+1)+0.02, col=rgb(0.98, 0.98, 0.98), border=rgb(0.98, 0.98, 0.98))
		}
		if (sample %in% samples.mlpa) source <- paste(c(source, "M"), collapse="+")
		
		axis(4, (numsamp-i+1.5), source, las=2, cex.axis=0.3, tck=-0.004, mgp=c(3, .15, 0), lwd.ticks=0.5)
	}
	
	# vertical grid lines separating chromosomes
	for (i in sizes$offset[2:(length(sizes$offset))]) {
		lines(c(i, i), c(0, numsamp)+1, type="l", col="lightgray", lty=1, lwd=0.7)	
	}
	lines(c(sum(sizes$size), sum(sizes$size)), c(0, numsamp)+1, type="l", col="lightgray", lty=1, lwd=1)
	
	# vertical grid lines indicating gene marks
	for (i in 1:nrow(markers)) {
		lines(c(markers$offset[i], markers$offset[i]), c(0, numsamp)+1, type="l", col="lightgray", lty=2, lwd=0.2)	
	}
	
	for (i in c(1:numsamp)) {
		sample <- samples[i]
		
		# get cnvs from array, WES, and MLPA
		if (!sample %in% c(samples.lowqual, unique(d.array$sample))) {
			cnvs <- d[d$sample.name == sample & d$source %in% c("WES", "MLPA"),]
		}
		else {
			cnvs <- d[d$sample.name == sample & d$source %in% c("array", "MLPA"),] # only plot non-WES calls (=SNParray & MLPA)		
		}
		
		# plot cnvs
		if (nrow(cnvs) > 0) {
			cnvs <- cnvs[order(cnvs$copy.count, cnvs$end-cnvs$start, decreasing=T),] # plot gains, larger cnvs first
			for (j in c(1:nrow(cnvs))) {
				width <- cnvs[j,"end.cum"]-cnvs[j,"start.cum"] 
				if (width<min.width) {
					cnvs[j,"start.cum"] = cnvs[j,"start.cum"] - (min.width-width)/2
					cnvs[j,"end.cum"] = cnvs[j,"end.cum"] + (min.width-width)/2
				}
				rect(cnvs[j,"start.cum"], (numsamp-i+1)+0.98, cnvs[j,"end.cum"], (numsamp-i+1)+0.02, col=cnvs$color[j], border=NA)	
			}			
		}	
		
	}
	
	# horizontal grid lines
	lastcase <- ""
	for (i in c(1:numsamp)) {
		case <- gsub("(D|R\\d?)$", "", rev(samples)[i])
		timepoint <- gsub("(.*)(D|R\\d?)$", "\\2", rev(samples)[i])
		if (case != lastcase) {
			lastcase <- case
			lines(c(0,sum(sizes$size)), c(i, i), type="l", col="lightgray", lty=1, lwd=1)	
			#axis(2, i+1, case, las=2, tick=F)
		}
		else {
			lines(c(0,sum(sizes$size)), c(i, i), type="l", col="lightgray", lty=2, lwd=0.5)	
		}
	}
	
	# axes
	axis(1,sizes$offset+sizes$size/2, sizes$chr, cex.axis=0.45, tick=F, las=3, mgp=c(3,.2,0), tck=-0.007, lwd.ticks=0.7, hadj=1)
	axis(3,markers$offset, markers$name, las=2, cex.axis=0.6, padj=markers$padj[order(markers$offset)], mgp=c(3, .3, 0), tck=-0.007, lwd.ticks=0.7)
	axis(2, (1:numsamp)+0.5, gsub("R\\d+$", "RR", gsub("(.*)(D|R\\d?)$", "\\1 \\2", rev(samples))), las=2, tick=F, cex.axis=0.6, mgp=c(3, .3, 0))
	box()
	legend(sum(sizes$size)/2, -1, c(0,1,3,4,"5+","n/a"), fill=c(col.0, col.1, col.3, col.4, col.5, "lightgray"), horiz=T, xjust=0.5, seg.len=10, xpd=T, cex=0.7)
}

# ---------------------------------------
# relapsing cases
# ---------------------------------------
pdf("~/p2ry8-crlf2/results/figures/cnv-plot-relapsing.pdf", paper="a4")
par(mar=c(4.5,6,6,3))
numsamp <- length(samples.matched)
plot.cnvs(samples.matched)
dev.off()


# ---------------------------------------
# non-relapsing cases (diagnosis only)
# ---------------------------------------
pdf("~/p2ry8-crlf2/results/figures/cnv-plot-nonrelapsing.pdf", paper="a4")
par(mar=c(4.5,6,6,3))
numsamp <- length(samples.diaonly)
plot.cnvs(samples.diaonly)
dev.off()
