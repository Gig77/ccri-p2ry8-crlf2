#library("RColorBrewer")
#X11.options(type="Xlib")

library(lattice)
library(gridExtra)
library(reshape)

# xenografts diagnosis: m1977-G-dia m1966-Y-dia m1035-108-dia m252-379-dia m1041-737-dia m247-833-dia m1059-92-dia m1037-839-dia m248-841-dia
# xenografts relapse: m1963-545-rel m1957-715-rel m1967-Y-rel m1060-108-rel m1069-737-rel
# 2nd xenograft relapse: m1964-545-rel 

# 737 hat im 2. relapse nur 8% blasten, deshalb ausgenommen von der klonalen analyse
#p <- "737" ; sample.RR <- "rem_rel3" ; exclude.chr <- c("X", "Y") ; genes.to.label <- c("CREBBP", "KRAS", "NRAS", "FLT3", "PTPN11", "FGFR2", "EPHA5") ; ; max.af.plot <- 0.7
#p <- "108" ; sample.RR <- "rem_rel2" ; exclude.chr <- NA ; genes.to.label <- c("JAK2", "KMT2D", "SI", "MEGF6", "NBEA", "RPS5", "BIRC6", "SLC26A2") ; max.af.plot <- 0.7
#p <- "AL9890" ; sample.RR <- "rem_rel2" ; exclude.chr <- NA ; genes.to.label <- c("CRYZ", "SEC16B", "FAT1") ; max.af.plot <- 0.7
#p <- "S23" ; sample.RR <- "rem_rel3" ; exclude.chr <- NA ; genes.to.label <- c("NPY2R", "CREBBP", "GPR63") ; max.af.plot <- 1

label.genes <- c("PAR1", "21", "X", "Y", "IKZF1", "IKZF2", "PAX5", "EBF1", "ETV6", "JAK2", "JAK3", "JAK1", "IL7R", "SYK", "CRLF2", "KRAS", "NRAS", "PTPN11", "CBL", "FLT3", "CDKN2A", "BTG1", "E2F3", "FOXP1", "AFF1", "XBP1", "TRRAP", "SETD2", "CREBBP", "KMT2D", "EP300", "USP9X", "BRCA2", "MET", "MSH6")

input <- list(
		"108" = list(
				title.suffix = "",
				timepoints =   c("D",           "DX",                    "R",           "RX",                    "RR"),
				identifiers =  c("108 rem_dia", "m1035-108-dia rem_xeno", "108 rem_rel", "m1060-108-rel rem_xeno", "108 rem_rel2"),
				xpos =         c(1,             2,                        4,             5,                        7),
				xenograft =    c(FALSE,         TRUE,                     FALSE,         TRUE,                     FALSE),
				transplant =   c(FALSE,         FALSE,                    FALSE,         FALSE,                    TRUE),
				blast.counts = c("91",          "NA",                     "93",          "NA",                     "85"),
				genes.to.label = label.genes,
				exclude.chr = c("MT"),
				max.af.plot = 1.0),
		"715" = list(
				title.suffix = "HD",
				timepoints =   c("D",           "R",           "RX",                    "RR"),
				identifiers =  c("715 rem_dia", "715 rem_rel", "m1957-715-rel rem_xeno", "715 rem_rel3"),
				xpos =         c(1,             3,             4,                        6),
				xenograft =    c(FALSE,         FALSE,         TRUE,                     FALSE),
				transplant =   c(FALSE,         FALSE,         FALSE,                    TRUE),
				blast.counts = c("92",          "92",          "NA",                     "92"),
				genes.to.label = label.genes,
				exclude.chr = c("MT"),
				max.af.plot = 1.0),
		"737" = list(
				title.suffix = "dic (9,20)",
				timepoints =   c("D",           "DX",                    "R",           "RX",                      "RR"),
				identifiers =  c("737 rem_dia", "m1041-737-dia rem_xeno", "737 rem_rel", "m1069-737-rel rem_xeno", "737 rem_rel3"),
				xpos =         c(1,             2,                        4,             5,                        7),
				xenograft =    c(FALSE,         TRUE,                     FALSE,         TRUE,                     FALSE),
				transplant =   c(FALSE,         FALSE,                    FALSE,         FALSE,                    TRUE),
				blast.counts = c("97",          "NA",                     "96",          "NA",                     "70"),
				genes.to.label = label.genes,
				exclude.chr = c("MT"),
				max.af.plot = 1.0),
#	"545" = list(
#   title.suffix = "HD",
#	  timepoints =   c("D",           "R",           "RX1",                   "RX2"),
#	  identifiers =  c("545 rem_dia", "545 rem_rel", "m1963-545-rel rem_xeno", "m1964-545-rel rem_xeno"),
#	  xenograft =    c(FALSE,         FALSE,         TRUE,                     TRUE),
#	  transplant =   c(FALSE,         FALSE,         FALSE,                    FALSE),
#	  blast.counts = c("97",           "86",         "NA",                      "NA"),
#	  genes.to.label = label.genes,
#	  exclude.chr = c("MT"),
#	  max.af.plot = 1),
		"839" = list(
				title.suffix = "iAMP21",
				timepoints =   c("D",           "DX",                    "R"),
				identifiers =  c("839 rem_dia", "m1037-839-dia rem_xeno", "839 rem_rel"),
				xpos =         c(1,             2,                        4),
				xenograft =    c(FALSE,         TRUE,                     FALSE),
				transplant =   c(FALSE,         FALSE,                    FALSE),
				blast.counts = c("92",           "NA",                    "86"),
				genes.to.label = label.genes,
				exclude.chr = c("MT"),
				max.af.plot = 1.0),
		"92" = list(
				title.suffix = "iAMP21",
				timepoints =   c("D",          "DX",                    "R"),
				identifiers =  c("92 rem_dia", "m1059-92-dia rem_xeno", "92 rem_rel"),
				xpos =         c(1,             2,                      4),
				xenograft =    c(FALSE,         TRUE,                   FALSE),
				transplant =   c(FALSE,         FALSE,                  FALSE),
				blast.counts = c("97",          "NA",                   "90"),
				genes.to.label = label.genes,
				exclude.chr = c("MT"),
				max.af.plot = 1.0),
		"841" = list(
				title.suffix = "HD",
				timepoints =   c("D",           "DX"),
				identifiers =  c("841 rem_dia", "m248-841-dia rem_xeno"),
				xpos =         c(1,             2),
				xenograft =    c(FALSE,         TRUE),
				transplant =   c(FALSE,         FALSE),
				blast.counts = c("98",          "NA"),
				genes.to.label = label.genes,
				exclude.chr = c("MT"),
				max.af.plot = 1.0),
		"379" = list(
				title.suffix = "",
				timepoints =   c("D",           "DX"),
				identifiers =  c("379 rem_dia", "m252-379-dia rem_xeno"),
				xpos =         c(1,             2),
				xenograft =    c(FALSE,         TRUE),
				transplant =   c(FALSE,         FALSE),
				blast.counts = c("97",          "NA"),
				genes.to.label = label.genes,
				exclude.chr = c("MT"),
				max.af.plot = 1.0),
		"833" = list(
				title.suffix = "",
				timepoints =   c("D",           "DX"),
				identifiers =  c("833 rem_dia", "m247-833-dia rem_xeno"),
				xpos =         c(1,             2),
				xenograft =    c(FALSE,         TRUE),
				transplant =   c(FALSE,         FALSE),
				blast.counts = c("96",          "NA"),
				genes.to.label = label.genes,
				exclude.chr = c("MT"),
				max.af.plot = 1.0),
		"AL9890" = list(
				title.suffix = "",
				timepoints =   c("D",             "R",               "RR"),
				identifiers =  c("AL9890 rem_dia", "AL9890 rem_rel", "AL9890 rem_rel2"),
				xpos =         c(1,                2,                3),
				xenograft =    c(FALSE,            FALSE,            FALSE),
				transplant =   c(FALSE,            FALSE,            FALSE),
				blast.counts = c("NA",             "NA",             "86"),
				genes.to.label = label.genes,
				exclude.chr = c("MT"),
				max.af.plot = 1.0),
		"S23" = list(
				title.suffix = "",
				timepoints =   c("D",           "R",           "RR"),
				identifiers =  c("S23 rem_dia", "S23 rem_rel", "S23 rem_rel3"),
				xpos =         c(1,              2,            3),
				xenograft =    c(FALSE,          FALSE,        FALSE),
				transplant =   c(FALSE,          FALSE,        FALSE),
				blast.counts = c("NA",           "NA",         "96"),
				genes.to.label = label.genes,
				exclude.chr = c("MT"),
				max.af.plot = 1.0)
)

min.cov <- 0
cov.max.std.dev <- 10
max.af <- 1
min.af <- 0

# mutations
mut.prim <- read.csv("/mnt/projects/p2ry8-crlf2/results/filtered-variants.cosmic.tsv", sep="\t")
mut.xeno <- read.csv("/mnt/projects/p2ry8-crlf2/results/filtered-variants.xenografts.cosmic.tsv", sep="\t")
mut <- rbind(mut.prim, mut.xeno)
mut <- mut[!grepl("^GL", mut$chr),]
mut <- mut[mut$status != "REJECT",]
mut <- mut[mut$deleterious=="yes",]
mut <- mut[mut$gene %in% label.genes,]
mut <- mut[,c("patient", "sample", "chr", "pos", "ref", "alt", "gene", "freq_leu", "cosmic_hits_aa")]
mut$freq_leu.raw <- mut$freq_leu

# adjust AF for blast counts
for (p in names(input)) {
	for (i in 1:length(input[[p]]$identifiers)) {
		sample <- gsub(paste0(p, " "), "", input[[p]]$identifiers[i])
		blasts <- input[[p]]$blast.counts[i]
		if (blasts != "NA") {
			mut$freq_leu[mut$patient == p & mut$sample == sample] <- mut$freq_leu.raw[mut$patient == p & mut$sample == sample] / as.numeric(blasts) * 100
		}
	}
}
mut$freq_leu[mut$freq_leu>1] <- 1

# cnvs
cnv <- read.csv("/mnt/projects/p2ry8-crlf2/data/cnvs_kinetics_figure.txt", sep="\t")

#mut <- mut[mut$dp_leu_tot >= min.cov & mut$dp_rem_tot >= min.cov,]  # minimum coverage filter
#mut <- mut[mut$dp_rem_tot < mean(mut$dp_rem_tot) + cov.max.std.dev * sd(mut$dp_rem_tot),] # maximum coverage filter

cols.merge <- c("chr", "pos", "ref", "alt", "gene", "cosmic_hits_aa")
col.af <- "freq_leu"

pdf(paste0("/mnt/projects/p2ry8-crlf2/results/figures/clonal-kinetics.pdf"), width=8, height=6)

for (p in names(input)) {
	
	mut.patient <- NULL
	num.tp <- length(input[[p]]$timepoints)
	x <- 1.8+1.5*(input[[p]]$xpos-1)
	
	# merge allelic frequencies from all timepoints
	
	for(i in 1:length(input[[p]]$timepoints)) {
		id <- input[[p]]$identifiers[i]
		tp <- input[[p]]$timepoints[i]
		mut.tp <- mut[paste(mut$patient, mut$sample) == id, c(cols.merge, col.af)]
		names(mut.tp)[names(mut.tp)==col.af] <- tp
		
		if (is.null(mut.patient)) {
			mut.patient <- mut.tp
		} else {
			mut.patient <- merge(mut.patient, mut.tp, by=cols.merge, all=TRUE)
		}
		
		print(paste(id, ":", nrow(mut.tp), "variants"))
	}
	
	for(tp in input[[p]]$timepoints) {
		mut.patient[is.na(mut.patient[,tp]),tp] <- 0
	}
	
	# remove mutations with allelic frequency < 0.1 in ALL samples
	mut.patient <- mut.patient[rowSums(mut.patient[,input[[p]]$timepoints] >= 0.1) > 0,]
	
	# keep only mutations seen in at least one primary sample or known hotspot mutations reported multiple times in COSMIC (excludes possible mouse polymorphisms)
	mut.patient <- mut.patient[apply(mut.patient[,input[[p]]$timepoints][!input[[p]]$xenograft], 1, max) > 0 | !mut.patient$cosmic_hits_aa %in% c("non-coding", "0", "1", "2", "3", "4"),]
	
	if (nrow(mut.patient) == 0) {
		plot(0, 0, main=paste0(p, ": No mutations found"), xlab="", ylab="", xaxt='n', yaxt='n')
		next
	}
	
	# color mutations by gene
	mut.genes <- as.factor(as.character(mut.patient$gene))
	colors <- c("#0080ff", "#0000ff", "darkgreen", "#ff0000", "orange", "#00ff00", "brown")
	mut.patient$color <- colors[as.integer(mut.genes) %% (length(colors)-1) + 1]
	
	# all mutations in one plot
	
	plot(0, 0, xlim=c(0.5, max(x)+0.6), ylim=c(-0.08-0.03*length(unique((cnv$gene[cnv$patient==p]))), input[[p]]$max.af.plot), type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="CNA                 Mutations (adj. allelic frequency)", main=ifelse(input[[p]]$title.suffix != "", paste(p, input[[p]]$title.suffix), p))
	abline(h=seq(0.1, 1, 0.1), lty=2, col="#eeeeee", lwd=1)
	abline(v=x, lty=2, col="black", lwd=1)
	abline(h=0.3, lty=1, col="#bbbbbb", lwd=2)
	abline(h=-0.05)
	box()
	axis(1, at = x, labels=input[[p]]$timepoints, padj=0.5)
	axis(2, at = seq(0, 1, 0.1), las = 1)
	legend("topleft", levels(mut.genes), cex = 0.8, fill = colors[(1:length(levels(mut.genes))) %% (length(colors)-1) + 1], bg = "white")
	
	# mutations
	
	for(i in 1:nrow(mut.patient)) {
		
		# primary kinetics
		if (max(mut.patient[i,input[[p]]$timepoints][!input[[p]]$xenograft]) > 0) {
			lines(x[!input[[p]]$xenograft], mut.patient[i,input[[p]]$timepoints][!input[[p]]$xenograft], type="b", col=mut.patient$color[i], pch=ifelse(mut.patient[i,input[[p]]$timepoints][!input[[p]]$xenograft] > 0, 19, 1), lwd=3.5, cex=1.5)
		}
		
		# xenograft kinetics
		for (s in which(input[[p]]$xenograft)) {
			xeno.col <- match(input[[p]]$timepoints[s], names(mut.patient))
			if (max(mut.patient[i,(xeno.col-1):xeno.col]) == 0) next
			lines(x[(s-1):s], mut.patient[i,(xeno.col-1):xeno.col], type="b", col=mut.patient$color[i], pch=ifelse(mut.patient[i,(xeno.col-1):xeno.col] > 0, 19, 1), lwd=1.5, lty=1, cex=1.5)
		}
		if (mut.patient$gene[i] != "" & mut.patient$gene[i] %in% input[[p]]$genes.to.label) {
			#text(1.5, mut.patient[i,input[[p]]$timepoints[1]], paste0(mut.patient$gene[i], " --"), adj=c(1, 0.5), cex=0.8)
			#text(max(x)+0.35, mut.patient[i,tail(input[[p]]$timepoints,1)], paste0("-- ", mut.patient$gene[i]), adj=c(0, 0.5), cex=0.8)
		}
	}	
	rect(par("usr")[1]+0.02, -0.048, par("usr")[2]-0.02, 0.098, col="#eeeeeeaa", border=NA)
	#lines(c(par("usr")[1]+0.02, par("usr")[2]-0.02), c(0.097, 0.097), col="#aaaaaa", lwd=0.5)
	for(i in 1:length(x)) {
		lines(c(x[i], x[i]), c(0.1, -0.05), lty=2, col="black", lwd=1)
	}
	
	# cnvs
	cnv.patient <- cnv[paste(cnv$patient, cnv$sample) %in% input[[p]]$identifiers,]
	y <- -0.1
	for (gene in unique(cnv.patient$gene)) {
		lines(c(x[1], x[length(x)]), c(y, y), lty=2, col="#cccccc", lwd=1)
		for (i in which(cnv.patient$gene==gene)) {
			print(paste(cnv.patient$patient[i], cnv.patient$sample[i]))
			print(input[[p]]$xpos[match(paste(cnv.patient$patient[i], cnv.patient$sample[i]), input[[p]]$identifiers)])
			xpos <- x[match(paste(cnv.patient$patient[i], cnv.patient$sample[i]), input[[p]]$identifiers)]
			cnv.col <- ifelse(cnv.patient$event[i] %in% c("deletion", "IK6"), "#0000ff", "red")
			if (cnv.patient$gene[i] == "PAR1") cnv.col <- "darkblue"
			rect(xpos-0.5, y+0.012, xpos+0.5, y-0.012, col=cnv.col, border=NA)
			text(1.2, y, gene, adj=c(1, 0.5), cex=0.65)
		}
		y <- y - 0.03
	}
	
	next
	
	# lattice plot with mutations split by gene
	
	m <- melt(data.patient, id.vars=c(cols.merge, "color"))
	m$mutation <- as.factor(paste0(m$chr, ":", m$pos, ":", m$ref, ">", m$alt))
	gene.colors <- m$color[match(levels(m$mutation), m$mutation)] # arrange colors by grouping factor
	
	print(xyplot(value ~ variable | gene, 
					group=mutation, 
					data=m, 
					type='b',
					par.settings=list(superpose.line = list(col=gene.colors), superpose.symbol = list(pch=19, col=gene.colors)),
					main=p, 
					as.table=TRUE, 
					par.strip.text=list(cex=0.7),
					xlab="Sample",
					ylab="Allelic frequency",
					ylim=c(-0.1,1.1),
					panel = function(...) {
						panel.abline(h=seq(0, 1, 0.1), col.line="#eeeeee", lty=2, lwd=0.5)
						panel.xyplot(...)
					}
			))
	
	#plot.conserved <- xyplot(value ~ variable | gene, group=mutation, data=m[m$group=="conserved",], type='b', layout=c(1, 4))
	#plot.unassigned <- xyplot(value ~ variable | gene, group=mutation, data=m[m$group=="unassigned",], type='b', layout=c(2, 4, 1))
	#grid.arrange(plot.conserved, plot.unassigned, ncol=length(unique(m$group)))
	
	#write.table(m, file=paste0("/mnt/projects/p2ry8-crlf2/results/figures/clonal-kinetics-", p, ".tsv"), col.names=T, row.names=F, sep="\t", quote=F)
}

dev.off()
