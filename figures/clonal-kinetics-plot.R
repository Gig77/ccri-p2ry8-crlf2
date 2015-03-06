#library("RColorBrewer")
#X11.options(type="Xlib")

#p <- "737" ; sample.RR <- "rem_rel3" ; exclude.chr <- c("X", "Y") ; genes.to.label <- c("CREBBP", "KRAS", "NRAS", "FLT3", "PTPN11", "FGFR2", "EPHA5")
p <- "108" ; sample.RR <- "rem_rel2" ; exclude.chr <- NA ; genes.to.label <- c("JAK2", "KMT2D", "BIRC6", "NBEA", "GNAQ", "GLRX", "CREBBP", "KRAS", "NRAS", "FLT3", "PTPN11", "MLL2")

min.cov <- 0
cov.max.std.dev <- 10
max.af <- 1
min.af <- 0
blast.count <- list("737.dia" = 97, "737.rel" = 96, "737.rel2" = 70, "108.dia" = 91, "108.rel" = 93, "108.rel2" = 85)

# exome seq variants
data <- read.csv("~/p2ry8-crlf2/results/filtered-variants.cosmic.tsv", sep="\t")
data <- data[data$patient == p & data$status != "REJECT",]
data <- data[data$non_silent==1,]  # only non-silent
#data <- data[data$var_type=="snp",]  # variant type filter
data <- data[!(data$chr %in% exclude.chr),]  # sex chromosome filter
#data <- data[data$dp_leu_tot >= min.cov & data$dp_rem_tot >= min.cov,]  # minimum coverage filter
#data <- data[data$dp_rem_tot < mean(data$dp_rem_tot) + cov.max.std.dev * sd(data$dp_rem_tot),] # maximum coverage filter

dia <- data[data$sample=="rem_dia", c("chr", "pos", "ref", "alt", "gene", "freq_leu")]
names(dia)[6] <- "dia"
rel <- data[data$sample=="rem_rel", c("chr", "pos", "ref", "alt", "gene", "freq_leu")]
names(rel)[6] <- "rel"
rel2 <- data[data$sample==sample.RR, c("chr", "pos", "ref", "alt", "gene", "freq_leu")]
names(rel2)[6] <- "rel2"

# merge both
m <- merge(dia, rel, all.x=T, all.y=T)
m <- merge(m, rel2, all.x=T, all.y=T)

#par(mfrow = c(4, 5), mar=c(2,3,2,1))

m[is.na(m$dia), "dia"] <- 0
m[is.na(m$rel), "rel"] <- 0
m[is.na(m$rel2), "rel2"] <- 0
m[m$dia>1, "dia"] <- 1
m[m$rel>1, "rel"] <- 1
m[m$rel2>1, "rel2"] <- 1

# remove variants with AF < 10% at all timepoints (likely false positives)
#m <- m[m$dia>=min.af | m$rel>=min.af | m$rel2>=min.af,]

# remove variants with AF > 70% (likely inaccurate measurement)
#m <- m[m$dia<=max.af & m$rel<=max.af & m$rel2<=max.af,]

# assign colors depending on AFs
m$col <- "#000000"
m$col[m$dia>0.25 & m$rel>0.1 & m$rel2>0.25] <- "#C6C7C8"
m$col[m$dia>0.25 & m$rel==0 & m$rel2==0] <- "#4A8ECC"
m$col[m$dia>0.25 & m$rel>0 & m$rel<=0.25 & m$rel2==0] <- "#231F20"
m$col[m$dia==0 & m$rel>0.25 & m$rel2 > 0.25] <- "#C9252B"
m$col[m$dia==0 & m$rel>0.25 & m$rel2 == 0] <- "#76C58F"
m$col[m$dia==0 & m$rel<=0.25 & m$rel2 == 0] <- "#548A2F"
m$col[m$dia==0 & m$rel<=0.25 & m$rel2 > 0.25] <- "#F58E7D"
m$col[m$dia==0 & m$rel==0 & m$rel2 > 0.25] <- "#BCAFD6"
m$col[m$dia==0 & m$rel==0 & m$rel2 <= 0.25] <- "#FFF57C"

pdf(paste0("~/p2ry8-crlf2/results/figures/clonal-kinetics-", p, ".pdf"), width=12)
plot(0, 0, xlim=c(1, 4.6), ylim=c(0, 0.7), type="n", xaxt="n", yaxt="n", xlab="", ylab="Allelic frequency", main=paste(p, " (n=", nrow(m), ")", sep=""))
axis(1, at=c(1.3, 2.8, 4.3), labels=c(paste0("D\n(", blast.count[[paste0(p, ".dia")]], "% blasts)"), paste0("R\n(", blast.count[[paste0(p, ".rel")]], "% blasts)"), paste0("RR\n(", blast.count[[paste0(p, ".rel2")]], "% blasts)")), padj=0.5)
axis(2, at = seq(0, 1, 0.1), las = 1) 

for(i in 1:nrow(m)) {
	fdia <- m[i,"dia"]
	frel <- m[i,"rel"]
	frel2 <- m[i,"rel2"]
	lines(c(1.3, 2.8, 4.3), c(fdia, frel, frel2), type="b", col=m[i,"col"], lwd=ifelse(m[i, "gene"] %in% genes.to.label, 4, 1))
	
	if (m[i, "gene"] %in% genes.to.label) {
		text(1.25, fdia, paste0(m[i, "gene"], " --"), adj=c(1, 0.5), cex=0.8)
		text(4.35, frel2, paste0("-- ", m[i, "gene"]), adj=c(0, 0.5), cex=0.8)
	}
}	

write.table(m, file=paste0("~/p2ry8-crlf2/results/figures/clonal-kinetics-", p, ".tsv"), col.names=T, row.names=F, sep="\t", quote=F)

dev.off()
