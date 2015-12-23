library(foreach)
library(doMC)
registerDoMC(20) # number of cores for parallel execution

rm(list=ls())

min.basequal <- 20
min.reads.rev <- 2
min.reads.fwd <- 2
hotspots <- read.delim("/mnt/projects/p2ry8-crlf2/scripts/signaling-hotspots.tsv", colClasses = c("character", "character", "numeric"))

cases.relapsing.matched <- c("108", "92", "460", "545", "564", "715", "737", "839", "B36", "BB16", "DL2", "GI8", "GI13", "HV57", "HV80", "LU3", "MA5", "N7", "S23", "SN18", "DS10898", "VS14645", "SE15285", "KE17247", "AL9890", "GL11356")
cases.relapsing.diaonly <- c("1060", "BJ17183")
cases.nonrelapsing <- c("242", "360", "365", "379", "400", "506", "769", "802", "833", "887", "841", "903", "948", "957", "961", "1066", "1089", "HW11537", "KT14158", "TL14516")
#cases.relapsing.matched <- c("564")
#cases.relapsing.diaonly <- c()
#cases.nonrelapsing <- c()

check_base <- function(p, cohort, sample, gene, chr, pos, ref, base, qual, mutations) {
	for (b in c("A", "T", "G", "C")) {
		#print(paste0("ALT", ": ", sum(base==tolower(b) | base==toupper(b)), " ", mean(qual[base==tolower(b) | base==toupper(b)])))
		#print(paste0("REF", ": ", sum(base=="," | base=="."), " ", mean(qual[base=="," | base=="."])))
		qual.forward.mean <- mean(qual[base==toupper(b)])
		qual.reverse.mean <- mean(qual[base==tolower(b)])
		qual.forward.max <- max(qual[base==toupper(b)])
		qual.reverse.max <- max(qual[base==tolower(b)])
		if (sum(base==tolower(b)) >= min.reads.rev && sum(base==toupper(b)) >= min.reads.fwd && qual.forward.mean >= min.basequal && qual.reverse.mean >= min.basequal)
		{
			num.alt <- sum(toupper(base)==b)
			tot <- length(base)
			freq <- num.alt/tot
			print(paste0(p, "\t", cohort, "\t", sample, "\t", gene, "\t", chr, ":", pos, "\t", ref, "->", b, "\tALT:", num.alt, "\tTOT:", tot, "\tFREQ:", freq, "\tQUAL FWD:", qual.forward.mean, "\tQUAL REV:", qual.reverse.mean))
			mutations[nrow(mutations)+1,] <- c(p, cohort, sample, gene, chr, pos, ref, b, num.alt, tot, freq, qual.forward.mean, qual.reverse.mean, qual.forward.max, qual.reverse.max)
		}	
	}
	
	return(mutations)
}

process_bam <- function(p, cohort, sample, bam, gene, chr, pos, mutations) {
	cmd <- paste0("~/tools/samtools-0.1.19/samtools view -bh -F 1024 ", bam, " ", chr, ":", pos, "-", pos+1, " | ~/tools/samtools-0.1.19/samtools mpileup -q 1 -f /mnt/projects/generic/data/broad/human_g1k_v37.fasta - 2>/dev/null | grep ", pos)
	#print(cmd)
	r <- system(cmd, intern=T)
	if (length(r) == 0) {
		print(paste("WARNING: NO COVERAGE:", p, cohort, sample, bam, gene, chr, pos))
		return(mutations)
	}
	#print(r)
	
	ref <- strsplit(r, "\t")[[1]][3]
	base <- strsplit(r, "\t")[[1]][5]
	base <- gsub("\\^.", "", base, perl=T)
	base <- gsub("\\$", "", base, perl=T)
	base <- strsplit(base, "")[[1]]
	#print(base)
	qual <- as.numeric(charToRaw(strsplit(r, "\t")[[1]][6]))-33
	#print(qual)
	
	mutations <- check_base(p, cohort, sample, gene, chr, pos, ref, base, qual, mutations)
	
	return(mutations)
}

# ------------------------------------------------------------------------------------

result <- foreach(i=1:nrow(hotspots), .verbose = FALSE) %dopar% {
	mutations <- data.frame(patient=character(0), cohort=character(0), sample=character(0), gene=character(0), chr=character(0), pos=numeric(0), ref=character(0), alt=character(0), alt.reads=numeric(0), tot.reads=numeric(0), frequency=numeric(0), qual.forward.mean=numeric(0), qual.reverse.mean=numeric(), qual.forward.max=numeric(0), qual.reverse.max=numeric(), stringsAsFactors=F)
	for (p in cases.relapsing.matched) {
		mutations <- process_bam(p, "relapsing", "diagnosis", paste0("/mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_", p, "D_realigned.bam"), hotspots[i, "gene"], hotspots[i, "chr"], hotspots[i, "pos"], mutations)
		mutations <- process_bam(p, "relapsing", "relapse", paste0("/mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_", p, "R_realigned.bam"), hotspots[i, "gene"], hotspots[i, "chr"], hotspots[i, "pos"], mutations)
		mutations <- process_bam(p, "relapsing", "remission", paste0("/mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_", p, "C_realigned.bam"), hotspots[i, "gene"], hotspots[i, "chr"], hotspots[i, "pos"], mutations)
	}
	for (p in cases.relapsing.diaonly) {
		mutations <- process_bam(p, "relapsing", "diagnosis", paste0("/mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_", p, "D_realigned.bam"), hotspots[i, "gene"], hotspots[i, "chr"], hotspots[i, "pos"], mutations)
		mutations <- process_bam(p, "relapsing", "remission", paste0("/mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_", p, "C_realigned.bam"), hotspots[i, "gene"], hotspots[i, "chr"], hotspots[i, "pos"], mutations)
	}
	for (p in cases.nonrelapsing) {
		mutations <- process_bam(p, "non-relapsing", "diagnosis", paste0("/mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_", p, "D_realigned.bam"), hotspots[i, "gene"], hotspots[i, "chr"], hotspots[i, "pos"], mutations)
		mutations <- process_bam(p, "non-relapsing", "remission", paste0("/mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_", p, "C_realigned.bam"), hotspots[i, "gene"], hotspots[i, "chr"], hotspots[i, "pos"], mutations)
	}
	mutations
}
result <- do.call("rbind", result)
write.table(result, file="/mnt/projects/p2ry8-crlf2/results/hotspot-mutations.tsv", col.names=T, row.names=F, sep="\t", quote=F)
