options(warn=1)

patients <- c("108", "92", "715", "737", "839", "B36", "BB16", "DL2", "GI8", "GI13", "HV57", "HV80", "N7", "S23", "DS10898", "VS14645", "SE15285", "AL9890", "GL11356")
title.suffix <- list("839" = "iAMP21", "B36" = "iAMP21", "92" = "iAMP21", "715" = "HD", "737" = "dic (9,20)")
label.genes <- c("P-C", "+21", "+X", "+Y", "-Y", "IKZF1", "IKZF2", "PAX5", "EBF1", "ETV6", "JAK2", "JAK3", "JAK1", "IL7R", "SYK", "CRLF2", "KRAS", "NRAS", "PTPN11", "CBL", "FLT3", "CDKN2A", "BTG1", "E2F3", "FOXP1", "AFF1", "XBP1", "TRRAP", "SETD2", "CREBBP", "KMT2D", "EP300", "USP9X", "BRCA2", "MET", "MSH6")

t <- read.csv("/mnt/projects/p2ry8-crlf2/results/filtered-variants.cosmic.tsv", sep="\t", stringsAsFactors=F)
t <- t[t$status != "REJECT" & t$deleterious == "yes" & t$freq_leu >= 0.1, c("patient", "sample", "chr", "pos", "gene", "ref", "alt", "freq_leu")]

# add gains and losses

# P2RY8-CRLF2 fusion
for (p in patients) {
  t <- rbind(t, setNames(data.frame(p, "rem_dia", "X", 1581466, "P-C", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
}
for (p in c("N7", "DL2", "108", "HV80", "DS10898", "GI8", "VS14645", "839", "BB16", "GL11356", "AL9890", "715", "SE15285")) {
  t <- rbind(t, setNames(data.frame(p, "rem_rel", "X", 1581466, "P-C", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
}

# chr21 somatic
t <- rbind(t, setNames(data.frame("108", "rem_dia", "21", 1, "+21", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("108", "rem_rel", "21", 1, "+21", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("737", "rem_dia", "21", 1, "+21", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("737", "rem_rel", "21", 1, "+21", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("92", "rem_dia", "21", 1, "+21", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("92", "rem_rel", "21", 1, "+21", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("GI8", "rem_dia", "21", 1, "+21", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("B36", "rem_dia", "21", 1, "+21", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("839", "rem_dia", "21", 1, "+21", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("839", "rem_rel", "21", 1, "+21", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("715", "rem_dia", "21", 1, "+21", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("715", "rem_rel", "21", 1, "+21", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("S23", "rem_rel", "21", 1, "+21", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# sex chromosomes
t <- rbind(t, setNames(data.frame("N7", "rem_rel", "X", 1, "+X", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("108", "rem_dia", "X", 1, "+X", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("108", "rem_rel", "X", 1, "+X", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("737", "rem_rel", "Y", 1, "-Y", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("HV80", "rem_dia", "X", 1, "+X", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("B36", "rem_dia", "X", 1, "+X", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("B36", "rem_rel", "X", 1, "+X", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("GL11356", "rem_dia", "Y", 1, "+Y", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("AL9890", "rem_rel", "X", 1, "+X", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("715", "rem_dia", "X", 1, "+X", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("SE15285", "rem_dia", "X", 1, "+X", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("SE15285", "rem_rel", "X", 1, "+X", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("S23", "rem_rel", "X", 1, "+X", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# IK6
t <- rbind(t, setNames(data.frame("N7", "rem_dia", "7", 50344378, "IKZF1", "wt", "IK6", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("N7", "rem_rel", "7", 50344378, "IKZF1", "wt", "IK6", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("DL2", "rem_dia", "7", 50344378, "IKZF1", "wt", "IK6", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("DL2", "rem_rel", "7", 50344378, "IKZF1", "wt", "IK6", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("108", "rem_dia", "7", 50344378, "IKZF1", "wt", "IK6", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("108", "rem_rel", "7", 50344378, "IKZF1", "wt", "IK6", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("HV57", "rem_dia", "7", 50344378, "IKZF1", "wt", "IK6", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("HV57", "rem_rel", "7", 50344378, "IKZF1", "wt", "IK6", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("737", "rem_dia", "7", 50344378, "IKZF1", "wt", "IK6", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("737", "rem_rel", "7", 50344378, "IKZF1", "wt", "IK6", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("GI8", "rem_rel", "7", 50344378, "IKZF1", "wt", "IK6", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# IKAROS deletion
t <- rbind(t, setNames(data.frame("HV80", "rem_dia", "7", 50344378, "IKZF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("HV80", "rem_rel", "7", 50344378, "IKZF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("DS10898", "rem_dia", "7", 50344378, "IKZF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("DS10898", "rem_rel", "7", 50344378, "IKZF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("92", "rem_dia", "7", 50344378, "IKZF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("92", "rem_rel", "7", 50344378, "IKZF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("B36", "rem_rel", "7", 50344378, "IKZF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("GI13", "rem_rel", "7", 50344378, "IKZF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# CDKN2A deletion
t <- rbind(t, setNames(data.frame("108", "rem_dia", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("108", "rem_rel", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("737", "rem_dia", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("737", "rem_rel", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("HV80", "rem_dia", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("HV80", "rem_rel", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("836", "rem_dia", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("836", "rem_rel", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("BB16", "rem_dia", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("BB16", "rem_rel", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("SE15285", "rem_dia", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("SE15285", "rem_rel", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("S23", "rem_dia", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("S23", "rem_rel", "9", 21967751, "CDKN2A", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# ETV6 deletion
t <- rbind(t, setNames(data.frame("N7", "rem_dia", "12", 11802788, "ETV6", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("N7", "rem_rel", "12", 11802788, "ETV6", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("GI13", "rem_rel", "12", 11802788, "ETV6", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# PAX5 deletion
t <- rbind(t, setNames(data.frame("737", "rem_dia", "9", 36838531, "PAX5", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("737", "rem_rel", "9", 36838531, "PAX5", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("HV80", "rem_dia", "9", 36838531, "PAX5", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("HV80", "rem_rel", "9", 36838531, "PAX5", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("HV80", "rem_rel", "9", 36838532, "PAX5", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("N7", "rem_rel", "9", 36838531, "PAX5", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("DL2", "rem_rel", "9", 36838531, "PAX5", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("VS14645", "rem_rel", "9", 36838531, "PAX5", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# EBF1 deletion
t <- rbind(t, setNames(data.frame("839", "rem_dia", "5", 158122923, "EBF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("839", "rem_rel", "5", 158122923, "EBF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("GL11356", "rem_dia", "5", 158122923, "EBF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("GL11356", "rem_rel", "5", 158122923, "EBF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# CBL deletion
t <- rbind(t, setNames(data.frame("BB16", "rem_dia", "11", 119076986, "CBL", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("BB16", "rem_rel", "11", 119076986, "CBL", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# SETD2 deletion
t <- rbind(t, setNames(data.frame("AL9890", "rem_dia", "3", 47057898, "SETD2", "wt", "deletion", 0.92, stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("GI8", "rem_rel", "3", 47057898, "SETD2", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# CREBBP deletion
t <- rbind(t, setNames(data.frame("AL9890", "rem_dia", "16", 3775056, "CREBBP", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# KMT2D deletion
t <- rbind(t, setNames(data.frame("GI8", "rem_rel", "12", 49412758, "KMT2D", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# MET deletion
t <- rbind(t, setNames(data.frame("92", "rem_rel", "7", 116312459, "MET", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# BTG1 deletion
t <- rbind(t, setNames(data.frame("DL2", "rem_dia", "12", 92534054, "BTG1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("DL2", "rem_rel", "12", 92534054, "BTG1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("GL11356", "rem_dia", "12", 92534054, "BTG1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("VS14645", "rem_rel", "12", 92534054, "BTG1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# E2F3 deletion
t <- rbind(t, setNames(data.frame("92", "rem_dia", "6", 20402137, "E2F3", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("92", "rem_rel", "6", 20402137, "E2F3", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# FOXP1 gain
t <- rbind(t, setNames(data.frame("BB16", "rem_dia", "3", 71003865, "FOXP1", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("BB16", "rem_rel", "3", 71003865, "FOXP1", "wt", "gain", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# AFF1 deletion
t <- rbind(t, setNames(data.frame("SE15285", "rem_dia", "4", 87856154, "AFF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("SE15285", "rem_rel", "4", 87856154, "AFF1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# USP9X deletion
t <- rbind(t, setNames(data.frame("HV57", "rem_dia", "X", 40944888, "USP9X", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("HV57", "rem_rel", "X", 40944888, "USP9X", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# XBP1 deletion
t <- rbind(t, setNames(data.frame("737", "rem_dia", "22", 29190548, "XBP1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("737", "rem_rel", "22", 29190548, "XBP1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))
t <- rbind(t, setNames(data.frame("108", "rem_rel", "22", 29190548, "XBP1", "wt", "deletion", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t)))

# JAK2 LOH
t <- rbind(t, setNames(data.frame("839", "rem_rel", "9", 1, "JAK2", "wt", "LOH", runif(1, 0.4, 0.6), stringsAsFactors=F), names(t))) # LOH complete chromosome arm

pdf("/mnt/projects/p2ry8-crlf2/results/figures/freq-scatter.pdf", width=12, paper='A4r')
par(mfrow=c(3,4), mar=c(2,1,2,1.5), oma=c(2.5,4,0,0))
for(p in patients) {
	freqd <- t[t$patient==p & t$sample=="rem_dia", !names(t) %in% c("patient", "sample")]
	colnames(freqd)[6] = "freq_dia"
	freqr <- t[t$patient==p & t$sample=="rem_rel", !names(t) %in% c("patient", "sample")]
	colnames(freqr)[6] = "freq_rel"
	merged <- merge(freqd, freqr, all.x=T, all.y =T)
	merged$freq_dia[is.na(merged$freq_dia)]	= 0
	merged$freq_rel[is.na(merged$freq_rel)]	= 0
	title <- p
	if (!is.null(title.suffix[[p]])) title <- paste(title, title.suffix[[p]])
	title <- paste0(title, " (n=", length(merged$freq_dia), ")")
	plot(0, 0, col="#00000033", xlim=c(0,1), ylim=c(0,1), main=title, pch='', xlab='', ylab='', bty='n')
	abline(h=seq(0, 1, 0.1), lty=2, col="#eeeeee", lwd=0.5)
	abline(v=seq(0, 1, 0.1), lty=2, col="#eeeeee", lwd=0.5)
	box()
	points(merged$freq_dia, merged$freq_rel, col="#00000033", pch=19, cex=0.7)
	merged.label <- merged[merged$gene %in% label.genes,]
	if (nrow(merged.label) > 0) {
	  merged.label$col <- "green"
	  merged.label$col[merged.label$alt %in% c("IK6", "deletion", "gain", "LOH")] <- rgb(141/255, 180/255, 226/255) # blueish
	  points(merged.label$freq_dia, merged.label$freq_rel, col=merged.label$col, pch=19, cex=0.7)
	  text(merged.label$freq_dia+0.015, merged.label$freq_rel+0.015, merged.label$gene, adj = c(0, 0), cex=0.7)
	}
	if (match(p, patients) %% 12 == 0 || match(p, patients) == length(patients)) {
	  mtext("AF diagnosis", side=1, outer=T, cex=1, line=1)
	  mtext("AF relapse", side=2, outer=T, cex=1, line=2)
	}
}
dev.off()
