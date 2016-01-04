library("RColorBrewer")
library("ggplot2")

genes.jak <- c("JAK2", "JAK3", "JAK1", "IL7R", "SYK", "CRLF2")
genes.ras <- c("KRAS", "NRAS", "PTPN11", "CBL", "FLT3")

# TABLE: filtered-variants.tsv
# read and filter input data
t <- read.delim("/mnt/projects/p2ry8-crlf2/results/filtered-variants.tsv", stringsAsFactors=F)
t <- t[t$status != "REJECT" & t$freq_leu >= 0.1 & t$non_silent==1,]
t <- t[!t$patient %in% c("MA5", "BJ14367", "LU3", "SN18", "460", "545", "564", "957"),]
t <- t[!t$sample %in% c("rem_rel2", "rem_rel3"),]
t <- t[,c("patient", "sample", "cohort", "chr", "pos", "ref", "alt", "gene", "freq_leu")]

# add mutations found only by hotspot caller but not MuTect
t <- rbind(t, setNames(data.frame("VS14645", "rem_rel", "relapsing", "1", 65310517, "C", "T", "JAK1", 0.116, stringsAsFactors=F), names(t)))

t$class <- NA
t$class[t$gene %in% genes.ras] <- "RTK/Ras"
t$class[t$gene %in% genes.jak] <- "JAK/STAT"

m <- t[!is.na(t$class),]
m$Gene <- factor(as.character(m$gene), c(genes.jak, genes.ras))
m$class <- factor(m$class, c("JAK/STAT", "RTK/Ras"))

# adjust AF for BLAST count (if available)
c <- read.delim("/mnt/projects/p2ry8-crlf2/results/clinical/Clinical data_P-C_v6.txt", stringsAsFactors=F)
c$blasts_dia[c$blasts_dia=="73 pB"] <- 73
c$blasts_dia <- as.numeric(c$blasts_dia)
names(c)[names(c)=="blasts_rel..1st."] <- "blasts_rel"
c$blasts_rel[c$blasts_rel=="isoliertes ZNS Rezidiv"] <- NA
c$blasts_rel <- as.numeric(as.character(c$blasts_rel))
m <- merge(m, c[,c("patient_id", "blasts_dia")], by.x = "patient", by.y="patient_id", all.x=T)
m <- merge(m, c[,c("patient_id", "blasts_rel")], by.x = "patient", by.y="patient_id", all.x=T)
m$blasts_rel[m$patient=="715" & m$sample=="rem_rel3"] <- 92
m$blasts_rel[m$patient=="108" & m$sample=="rem_rel2"] <- 85
m$freq_leu.adj <- NA
m$freq_leu.adj[m$sample=="rem_dia"] <- with(m[m$sample=="rem_dia",], ifelse(!is.na(blasts_dia) & blasts_dia > 0, freq_leu / blasts_dia * 100, freq_leu))
m$freq_leu.adj[m$sample %in% c("rem_rel", "rem_rel2", "rem_rel3")] <- with(m[m$sample %in% c("rem_rel", "rem_rel2", "rem_rel3"),], ifelse(!is.na(blasts_rel) & blasts_rel > 0, freq_leu / blasts_rel * 100, freq_leu))
m$freq_leu.adj[is.na(m$freq_leu.adj)] <- m$freq_leu[is.na(m$freq_leu.adj)] 
m$freq_leu.adj[m$freq_leu.adj>=1] <- 0.99
m$freq_leu.adj[m$freq_leu.adj>0.9] <- m$freq_leu.adj[m$freq_leu.adj>0.9] / 2

m$sample <- as.character(m$sample)
m$sample[m$sample %in% c("rem_dia")] <- "Diagnosis" 
m$sample[m$sample %in% c("rem_rel", "rem_rel2", "rem_rel3")] <- "Relapse"
m$sample <- as.factor(m$sample)

# AF of only conserved mutations at diagnosis (will be marked with dots in plot)
m$freq_leu.adj.conserved <- ifelse((m$sample=="Diagnosis" & paste(m$patient, m$chr, m$pos, m$ref, m$alt) %in% with(m[m$sample=="Relapse",], paste(patient, chr, pos, ref, alt))) |
                                   (m$sample=="Relapse" & paste(m$patient, m$chr, m$pos, m$ref, m$alt) %in% with(m[m$sample=="Diagnosis",], paste(patient, chr, pos, ref, alt))), 
                                   m$freq_leu.adj, NA)

# assign colors to gene/cohort combinations
genes.col <- c(rep(brewer.pal(length(genes.jak)+3, "PuBu")[4:(length(genes.jak)+3)], 2), 
			   rep(brewer.pal(length(genes.ras)+3, "OrRd")[4:(length(genes.ras)+3)], 2))
names(genes.col) <- c(paste(genes.jak, "relapsing"), paste(genes.jak, "non-relapsing"),
					  paste(genes.ras, "relapsing"), paste(genes.ras, "non-relapsing"))

# assign shapes to gene/cohort combinations
genes.shape <- rep(NA, length(genes.col))
names(genes.shape) <- names(genes.col)
genes.shape[paste(genes.ras, "relapsing")] <- seq(0, length(genes.ras)-1) %% 3 + 15  # filled  
genes.shape[paste(genes.ras, "non-relapsing")] <- seq(0, length(genes.ras)-1) %% 3   # empty
genes.shape[paste(genes.jak, "relapsing")] <- seq(0, length(genes.jak)-1) %% 3 + 15  # filled  
genes.shape[paste(genes.jak, "non-relapsing")] <- seq(0, length(genes.jak)-1) %% 3   # empty

# combined factor gene + cohort
m$gene_cohort <- factor(paste(m$gene, m$cohort), levels=names(genes.col))

# manually compute horizontal value including jitter b/c we use the same jitter twice
set.seed(471)
m$x <- as.numeric(m$class) + runif(nrow(m),-0.2,0.2)

#m$shape <- as.integer(m$Gene) %% 3
#m$shape[m$cohort == "relapsing"] <- m$shape[m$cohort == "relapsing"] + 15

pdf("/mnt/projects/p2ry8-crlf2/results/figures/af-plot.pdf", width=10)
#pdf("/mnt/projects/p2ry8-crlf2/results/figures/af-plot-without-box.pdf", width=10)
p <- ggplot(m, aes(class, freq_leu.adj))
p <- p + theme_bw()
p <- p + facet_grid(.~sample, scales="free")
p <- p + geom_boxplot(colour=class, outlier.shape = NA)
p <- p + geom_point(aes(x=x, shape = gene_cohort, color = gene_cohort), size=3)
p <- p + geom_point(aes(x=x, y=freq_leu.adj.conserved), color="white", size=1)
#p <- p + scale_shape_identity()
#p <- p + geom_point(data=m[m$cohort=="relapsing",], aes(shape = Gene, color = Gene), position = position_jitter(width=0.2), size=3)
#p <- p + geom_point(data=m[m$cohort=="non-relapsing",], aes(shape = Gene, color = Gene), position = position_jitter(width=0.2), size=3)
p <- p + scale_shape_manual(values=genes.shape)
p <- p + scale_colour_manual(values=genes.col)
p <- p + scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0.1, 0.7))
p <- p + xlab("\nGroup") + ylab("Adj. allelic frequency\n")
print(p)
dev.off()

# print p-values
p.jakstat <- wilcox.test(m$freq_leu.adj[m$sample=="Diagnosis" & m$class=="JAK/STAT"], m$freq_leu.adj[m$sample=="Relapse" & m$class=="JAK/STAT"])$p.value
sprintf("JAK/STAT p-value adj. AF diagnosis vs. relapse: %.2g", p.jakstat)
p.ras <- wilcox.test(m$freq_leu.adj[m$sample=="Diagnosis" & m$class=="RTK/Ras"], m$freq_leu.adj[m$sample=="Relapse" & m$class=="RTK/Ras"])$p.value
sprintf("RTK/Ras p-value adj. AF diagnosis vs. relapse: %.2g", p.ras)


