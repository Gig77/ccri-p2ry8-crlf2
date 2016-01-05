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

label.genes <- c("P-C", "+21", "+X", "+Y", "-Y", "IKZF1", "IKZF2", "PAX5", "EBF1", "ETV6", "JAK2", "JAK3", "JAK1", "IL7R", "SYK", "CRLF2", "KRAS", "NRAS", "PTPN11", "CBL", "FLT3", "CDKN2A", "BTG1", "E2F3", "FOXP1", "AFF1", "XBP1", "TRRAP", "SETD2", "CREBBP", "KMT2D", "EP300", "USP9X", "BRCA2", "MET", "MSH6")

input <- list(
	"108" = list(
		timepoints =   c("D",           "Dxg",                    "R",           "Rxg",                    "RR"),
		identifiers =  c("108 rem_dia", "m1035-108-dia rem_xeno", "108 rem_rel", "m1060-108-rel rem_xeno", "108 rem_rel2"),
		xenograft =    c(FALSE,         TRUE,                     FALSE,         TRUE,                     FALSE),
		transplant =   c(FALSE,         FALSE,                    FALSE,         FALSE,                    TRUE),
		blast.counts = c("91",          "NA",                     "93",          "NA",                     "85"),
		genes.to.label = label.genes,
		exclude.chr = c("MT"),
		max.af.plot = 1),
	"715" = list(
	  timepoints =   c("D",           "R",           "Rxg",                    "RR"),
	  identifiers =  c("715 rem_dia", "715 rem_rel", "m1957-715-rel rem_xeno", "715 rem_rel3"),
	  xenograft =    c(FALSE,         FALSE,         TRUE,                     FALSE),
	  transplant =   c(FALSE,         FALSE,         FALSE,                    TRUE),
	  blast.counts = c("92",          "92",          "NA",                     "NA"),
	  genes.to.label = label.genes,
	  exclude.chr = c("MT"),
	  max.af.plot = 1),
	"737" = list(
    timepoints =   c("D",           "Dxg",                    "R",           "Rxg",                    "RR"),
    identifiers =  c("737 rem_dia", "m1041-737-dia rem_xeno", "737 rem_rel", "m1069-737-rel rem_xeno", "737 rem_rel3"),
    xenograft =    c(FALSE,         TRUE,                     FALSE,         TRUE,                     FALSE),
    transplant =   c(FALSE,         FALSE,                    FALSE,         FALSE,                    TRUE),
    blast.counts = c("97",          "NA",                     "96",          "NA",                     "NA"),
    genes.to.label = label.genes,
    exclude.chr = c("MT"),
    max.af.plot = 1),
	"545" = list(
	  timepoints =   c("D",           "R",           "Rxg1",                   "Rxg2"),
	  identifiers =  c("545 rem_dia", "545 rem_rel", "m1963-545-rel rem_xeno", "m1964-545-rel rem_xeno"),
	  xenograft =    c(FALSE,         FALSE,         TRUE,                     TRUE),
	  transplant =   c(FALSE,         FALSE,         FALSE,                    FALSE),
	  blast.counts = c("97",           "86",         "NA",                      "NA"),
	  genes.to.label = label.genes,
	  exclude.chr = c("MT"),
	  max.af.plot = 1),
	"839" = list(
	  timepoints =   c("D",           "Dxg",                    "R"),
	  identifiers =  c("839 rem_dia", "m1037-839-dia rem_xeno", "839 rem_rel"),
	  xenograft =    c(FALSE,         TRUE,                     FALSE),
	  transplant =   c(FALSE,         FALSE,                    FALSE),
	  blast.counts = c("92",           "NA",                    "86"),
	  genes.to.label = label.genes,
	  exclude.chr = c("MT"),
	  max.af.plot = 1),
	"92" = list(
	  timepoints =   c("D",          "Dxg",                   "R"),
	  identifiers =  c("92 rem_dia", "m1059-92-dia rem_xeno", "92 rem_rel"),
	  xenograft =    c(FALSE,         TRUE,                   FALSE),
	  transplant =   c(FALSE,         FALSE,                  FALSE),
	  blast.counts = c("97",          "NA",                   "90"),
	  genes.to.label = label.genes,
	  exclude.chr = c("MT"),
	  max.af.plot = 1),
	"841" = list(
	  timepoints =   c("D",           "Dxg"),
	  identifiers =  c("841 rem_dia", "m248-841-dia rem_xeno"),
	  xenograft =    c(FALSE,         TRUE),
	  transplant =   c(FALSE,         FALSE),
	  blast.counts = c("98",          "NA"),
	  genes.to.label = label.genes,
	  exclude.chr = c("MT"),
	  max.af.plot = 1),
	"379" = list(
	  timepoints =   c("D",           "Dxg"),
	  identifiers =  c("379 rem_dia", "m252-379-dia rem_xeno"),
	  xenograft =    c(FALSE,         TRUE),
	  transplant =   c(FALSE,         FALSE),
	  blast.counts = c("97",          "NA"),
	  genes.to.label = label.genes,
	  exclude.chr = c("MT"),
	  max.af.plot = 1),
	"833" = list(
	  timepoints =   c("D",           "Dxg"),
	  identifiers =  c("833 rem_dia", "m247-833-dia rem_xeno"),
	  xenograft =    c(FALSE,         TRUE),
	  transplant =   c(FALSE,         FALSE),
	  blast.counts = c("96",          "NA"),
	  genes.to.label = label.genes,
	  exclude.chr = c("MT"),
	  max.af.plot = 1)
)

min.cov <- 0
cov.max.std.dev <- 10
max.af <- 1
min.af <- 0

# read data
data.prim <- read.csv("/mnt/projects/p2ry8-crlf2/results/filtered-variants.cosmic.tsv", sep="\t")
data.xeno <- read.csv("/mnt/projects/p2ry8-crlf2/results/filtered-variants.xenografts.cosmic.tsv", sep="\t")
data <- rbind(data.prim, data.xeno)
data <- data[!grepl("^GL", data$chr),]
data <- data[data$status != "REJECT",]
data <- data[data$deleterious=="yes",]

#data <- data[data$dp_leu_tot >= min.cov & data$dp_rem_tot >= min.cov,]  # minimum coverage filter
#data <- data[data$dp_rem_tot < mean(data$dp_rem_tot) + cov.max.std.dev * sd(data$dp_rem_tot),] # maximum coverage filter

cols.merge <- c("chr", "pos", "ref", "alt", "gene", "non_silent", "deleterious", "cosmic_hits_nt", "dp_rem_tot")
col.af <- "freq_leu"

pdf(paste0("/mnt/projects/p2ry8-crlf2/results/figures/clonal-kinetics-xenografts.pdf"), width=12)

for (p in names(input)) {

  data.patient <- NULL
  num.tp <- length(input[[p]]$timepoints)
  x <- seq(1.3, 1.5*(num.tp-1)+1.3, by=1.5)
  
  # merge allelic frequencies from all timepoints
  
  for(i in 1:length(input[[p]]$timepoints)) {
    id <- input[[p]]$identifiers[i]
    tp <- input[[p]]$timepoints[i]
    data.tp <- data[paste(data$patient, data$sample) == id, c(cols.merge, col.af)]
    names(data.tp)[names(data.tp)==col.af] <- tp
    
    if (is.null(data.patient)) {
      data.patient <- data.tp
    } else {
      data.patient <- merge(data.patient, data.tp, by=cols.merge, all=TRUE)
    }
    
    print(paste(id, ":", nrow(data.tp), "variants"))
  }
  
  for(tp in input[[p]]$timepoints) {
    data.patient[is.na(data.patient[,tp]),tp] <- 0
  }

  # exclude mutations from unwanted chromosomes
  
  if (!is.null(input[[p]]$exclude.chr)) {
    data.patient <- data.patient[!(data.patient$chr %in% input[[p]]$exclude.chr),]  # sex chromosome filter to not distort AF extimates
  }
  
  # filter mutations based on conservation pattern
  # keep variants that are
  # a) conserved between at least two samples (counting xenografts as one) OR
  # b) predicted deleterious COSMIC variant with AF > 10% (excluding xenografts)
  
  primary.col <- input[[p]]$timepoints[!input[[p]]$xenograft]
  xeno.col <- input[[p]]$timepoints[input[[p]]$xenograft]
  primary <- rowSums(data.patient[, primary.col, drop=FALSE] > 0)
  xeno <- rowSums(data.patient[, xeno.col, drop=FALSE] > 0) 
  conserved <- primary > 1 | (primary > 0 & xeno > 0)
  cosmic <- data.patient$cosmic_hits_nt > 0 & data.patient$deleterious=="yes" & apply(data.patient[,input[[p]]$timepoints], 1, max) >= 0.1
  interesting <- data.patient$gene %in% label.genes & apply(data.patient[,input[[p]]$timepoints], 1, max) >= 0.1
  data.patient <- data.patient[(primary & interesting) | conserved,]
  
  if (nrow(data.patient) == 0) {
    plot(0, 0, main=paste0(p, ": No mutations found"), xlab="", ylab="", xaxt='n', yaxt='n')
    next
  }

  # color mutations by gene
  data.patient$color <- c("#0080ff", "#ff00ff", "darkgreen", "#ff0000", "orange", "#00ff00", "brown")[as.integer(as.factor(as.character(data.patient$gene))) %% 6 + 1]
  
  # define groups based on AF at each timepoint
  #group.conserved <- data.patient[rowSums(apply(data.patient[,input[[p]]$timepoints], 2, function(x) x > 0.2))>=4,]
  #group.notdia <- data.patient[data.patient$D == 0 & data.patient$Dxg == 0 & data.patient$R > 0.2 & data.patient$Rxg > 0.2 & data.patient$RR > 0.2,]
  data.patient$group <- "unassigned"
  
  # all mutations in one plot
  
  plot(0, 0, xlim=c(1, 1.5*(num.tp-1)+1.6), ylim=c(0, input[[p]]$max.af.plot), type="n", xaxt="n", yaxt="n", xlab="", ylab="Allelic frequency", main=paste(p, " (n=", nrow(data.patient), ")", sep=""))
  abline(h=seq(0, 1, 0.1), lty=2, col="#eeeeee", lwd=0.5)
  axis(1, at = x, labels=paste0(input[[p]]$timepoints, "\n(", input[[p]]$blast.counts, "% blasts)"), padj=0.5)
  axis(2, at = seq(0, 1, 0.1), las = 1) 
  
  # draw group polygons
  #polygon(c(x, rev(x)), c(apply(group.conserved[,input[[p]]$timepoints], 2, min), rev(apply(group.conserved[,input[[p]]$timepoints], 2, max))), col=rgb(1, 1, 0, 0.15), border=NA)
  #polygon(c(x, rev(x)), c(apply(group.notdia[,input[[p]]$timepoints], 2, min), rev(apply(group.notdia[,input[[p]]$timepoints], 2, max))), col=rgb(0, 0, 1, 0.15), border=NA)
  
  for(i in 1:nrow(data.patient)) {
    lines(x, data.patient[i,input[[p]]$timepoints], type="b", col=data.patient$color[i], pch=19, lwd=ifelse(data.patient$gene[i] %in% label.genes, 5, 1.5))
    if (data.patient$gene[i] != "" & data.patient$gene[i] %in% input[[p]]$genes.to.label) {
      text(1.25, data.patient[i,input[[p]]$timepoints[1]], paste0(data.patient$gene[i], " --"), adj=c(1, 0.5), cex=0.8)
      text(1.5*(num.tp-1)+1.35, data.patient[i,tail(input[[p]]$timepoints,1)], paste0("-- ", data.patient$gene[i]), adj=c(0, 0.5), cex=0.8)
    }
  }	
  
  # lattice plot with mutations split by gene

  m <- melt(data.patient, id.vars=c(cols.merge, "group", "color"))
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


