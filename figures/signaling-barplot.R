library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)

rm(list=ls())

cases.relapsing.matched <- c("108", "92", "460", "545", "564", "715", "737", "839", "B36", "BB16", "DL2", "GI8", "GI13", "HV57", "HV80", "LU3", "MA5", "N7", "S23", "SN18", "DS10898", "VS14645", "SE15285", "KE17247", "AL9890", "GL11356")
cases.relapsing.diaonly <- c("1060", "BJ17183")
cases.nonrelapsing <- c("242", "360", "365", "379", "400", "506", "769", "802", "833", "887", "841", "903", "948", "957", "961", "1066", "1089", "HW11537", "KT14158", "TL14516")
cases.exclude <- c("MA5", "BJ14367", "LU3", "SN18", "460", "545", "564", "957")

cols <- c(		"KRAS" = brewer.pal(8, "OrRd")[2],
				"KRAS:12:25378561:G>A" = brewer.pal(8, "OrRd")[2],
				"KRAS:12:25378562:C>T" = brewer.pal(8, "OrRd")[2],
				"KRAS:12:25378647:T>A" = brewer.pal(8, "OrRd")[2],
				"KRAS:12:25398281:C>T" = brewer.pal(8, "OrRd")[2],
				"KRAS:12:25398284:C>T" = brewer.pal(8, "OrRd")[2],
				"KRAS:12:25398284:C>A" = brewer.pal(8, "OrRd")[2],
				"KRAS:12:25398285:C>G" = brewer.pal(8, "OrRd")[2],
				"KRAS:12:25398285:C>T" = brewer.pal(8, "OrRd")[2],

				"NRAS" = brewer.pal(8, "OrRd")[4],
				"NRAS:1:115256529:T>C" = brewer.pal(8, "OrRd")[4],
				"NRAS:1:115258744:C>T" = brewer.pal(8, "OrRd")[4],
				"NRAS:1:115258747:C>T" = brewer.pal(8, "OrRd")[4],
				"NRAS:1:115258747:C>G" = brewer.pal(8, "OrRd")[4],
				"NRAS:1:115258748:C>T" = brewer.pal(8, "OrRd")[4],
		
				"PTPN11" = brewer.pal(8, "OrRd")[6],
				"PTPN11:12:112888210:G>A" = brewer.pal(8, "OrRd")[6],
				"PTPN11:12:112888211:A>T" = brewer.pal(8, "OrRd")[6],

				"FLT3" = brewer.pal(8, "OrRd")[8],
				"FLT3:13:28599081:C>A" = brewer.pal(8, "OrRd")[8],
				"FLT3:13:28592629:T>C" = brewer.pal(8, "OrRd")[8],
				
				"JAK1" = brewer.pal(6, "PuBu")[2],
				"JAK1:1:65305426:G>C" = brewer.pal(6, "PuBu")[2],
				"JAK1:1:65310517:C>T" = brewer.pal(6, "PuBu")[2],
				"JAK1:1:65344759:C>T" = brewer.pal(6, "PuBu")[2],
				
				"JAK2" = brewer.pal(6, "PuBu")[4],
				"JAK2:9:5073749:G>A" = brewer.pal(6, "PuBu")[4],
				"JAK2:9:5073753:T>C" = brewer.pal(6, "PuBu")[4],
				"JAK2:9:5078357:A>T" = brewer.pal(6, "PuBu")[4],
				"JAK2:9:5078360:A>G" = brewer.pal(6, "PuBu")[4],
				"JAK2:9:5078361:G>C" = brewer.pal(6, "PuBu")[4],
				"JAK2:9:5078362:A>T" = brewer.pal(6, "PuBu")[4],
				"JAK2:9:5078362:A>C" = brewer.pal(6, "PuBu")[4],
				"JAK2:9:5078395:C>G" = brewer.pal(6, "PuBu")[4],
				"JAK2:9:5090497:G>A" = brewer.pal(6, "PuBu")[4],

				"JAK3" = brewer.pal(6, "PuBu")[6],
				"JAK3:19:17945969:C>T" = brewer.pal(6, "PuBu")[6],
				"JAK3:19:17945970:G>A" = brewer.pal(6, "PuBu")[6],
				
				"CRLF2" = "black",
				"CRLF2:X:1314966:A>C" = "black",
				
				"diagnosis" = brewer.pal(8, "OrRd")[6],
				"relapse" = brewer.pal(6, "PuBu")[6])
				
# get KRAS, NRAS, PTPN11 mutations from hotspot mutation caller
m <- read.delim("/mnt/projects/p2ry8-crlf2/results/hotspot-mutations.tsv")

# add in MuTect calls in hotspot regions plus FLT3
hotspots <- read.delim("/mnt/projects/p2ry8-crlf2/scripts/signaling-hotspots.tsv", colClasses = c("character", "character", "numeric"))
mutect <- read.delim("/mnt/projects/p2ry8-crlf2/results/filtered-variants.cosmic.tsv", stringsAsFactors=F)
mutect <- mutect[(paste(mutect$chr, mutect$pos) %in% paste(hotspots$chr, hotspots$pos) | mutect$gene=="FLT3") & mutect$status != "REJECT" & mutect$non_silent==1, c("patient", "sample", "cohort", "chr", "pos", "ref", "alt", "gene", "dp_leu_tot", "dp_leu_var", "freq_leu")]
mutect$sample[mutect$sample=="rem_dia"] <- "diagnosis"
mutect$sample[mutect$sample=="rem_rel"] <- "relapse"
m <- merge(m, mutect, by=c("patient", "sample", "cohort", "gene", "chr", "pos", "ref", "alt"), all=T)
m$source <- NA
m$source[!is.na(m$frequency) & !is.na(m$freq_leu)] <- "hotspot+mutect"
m$source[is.na(m$frequency) & !is.na(m$freq_leu)] <- "mutect"
m$source[!is.na(m$frequency) & is.na(m$freq_leu)] <- "hotspot"
m$frequency[is.na(m$frequency)] <- m$freq_leu[is.na(m$frequency)]
m$alt.reads[is.na(m$alt.reads)] <- m$dp_leu_var[is.na(m$alt.reads)]
m$tot.reads[is.na(m$tot.reads)] <- m$dp_leu_tot[is.na(m$tot.reads)]
m <- m[,c("source", "patient", "cohort", "sample", "gene", "chr", "pos", "ref", "alt", "frequency")]

# exclude cases with subclonal P/C fusion at diagnosis
m <- m[!m$patient %in% cases.exclude,]

# add non-existent mutations with AF 0 to get complete legend in ggplot
m <- rbind(m, data.frame(source=NA, patient="108", cohort="relapsing", sample="diagnosis", gene="JAK1", chr=1, pos=65305426, ref="G", alt="C", frequency=0))

# add columns
m$mut <- paste0(m$gene, ":", m$chr, ":", m$pos, ":", m$ref, ">", m$alt)
m$mut.short <- paste0(m$ref, ">", m$alt)

# normalize by blast counts
c <- read.delim("/mnt/projects/p2ry8-crlf2/results/clinical/Clinical data_P-C_v6.txt", stringsAsFactors=F)
names(c)[names(c)=="blasts_dia"] <- "blasts.dia"
names(c)[names(c)=="blasts_rel..1st."] <- "blasts.rel"
c$blasts.dia[c$blasts.dia=="73 pB"] <- 73
c$blasts.dia <- suppressWarnings(as.numeric(c$blasts.dia))
c$blasts.rel <- suppressWarnings(as.numeric(as.character(c$blasts.rel)))
m <- merge(m, c[,c("patient_id", "blasts.dia")], by.x = "patient", by.y="patient_id", all.x=T)
m <- merge(m, c[,c("patient_id", "blasts.rel")], by.x = "patient", by.y="patient_id", all.x=T)
m$frequency.norm <- NA
m$frequency.norm[m$sample=="diagnosis"] <- with(m[m$sample=="diagnosis",], ifelse(!is.na(blasts.dia) & blasts.dia > 0, frequency / blasts.dia * 100, frequency))
m$frequency.norm[m$sample=="relapse"] <- with(m[m$sample=="relapse",], ifelse(!is.na(blasts.rel) & blasts.rel > 0, frequency / blasts.rel * 100, frequency))
m$frequency.norm[is.na(m$frequency.norm)] <- m$frequency[is.na(m$frequency.norm)] 
m$frequency.norm[m$frequency.norm>=1] <- 0.99

# label cases from matched cohort
m$patient.label <- as.factor(ifelse(m$patient %in% cases.relapsing.matched, paste0(m$patient,"*"), as.character(m$patient)))

# subset samples
m.relapsing <- m[m$cohort=="relapsing",]
write.table(m[!is.na(m$source),], "/mnt/projects/p2ry8-crlf2/results/figures/signaling-barplot.mutations.tsv", row.names=F, col.names=T, sep="\t", quote=F)

m.relapsing.dia <- m.relapsing[m.relapsing$sample=="diagnosis",]
m.relapsing.dia <- ddply(m.relapsing.dia[order(m.relapsing.dia$frequency.norm, decreasing=T),], .(patient), transform, midpoint=cumsum(frequency.norm)-0.5*frequency.norm+0.005) # add midpoint for label plotting
m.relapsing.rel <- m.relapsing[m.relapsing$sample=="relapse",]
m.relapsing.rel <- ddply(m.relapsing.rel[order(m.relapsing.rel$frequency.norm, decreasing=T),], .(patient), transform, midpoint=cumsum(frequency.norm)-0.5*frequency.norm+0.005) # add midpoint for label plotting
m.relapsing.cons <- merge(m.relapsing.dia, m.relapsing.rel, by=c("patient", "mut"))[,c("patient", "mut", "frequency.norm.x", "frequency.norm.y")]
m.relapsing.cons <- paste0(m.relapsing.cons$patient, ":", m.relapsing.cons$mut)
m.nonrelapsing <- m[m$cohort=="non-relapsing" & m$sample=="diagnosis",]

pdf("/mnt/projects/p2ry8-crlf2/results/figures/signaling-barplot.relapsing.pdf", width=18, height=13)

#---
# RELAPSING, DIA
#---

# reorder patients by cumulative frequency
sorted <- aggregate(frequency.norm~patient.label, data=m.relapsing.dia, FUN=sum)
sorted <- sorted[order(sorted$frequency.norm, decreasing=T), c("patient.label", "frequency.norm")]
m.relapsing.dia$patient.label <- factor(m.relapsing.dia$patient.label, sorted$patient.label)

# find heterogenious samples
het.dia <- aggregate(frequency.norm~patient, data=m.relapsing.dia, FUN=length)
het.dia.patient <- het.dia[het.dia$frequency.norm>1,"patient"]
#m.relapsing.dia$group <- ifelse(m.relapsing.dia$patient %in% het.dia.patient, sprintf("heterogenious (%.0f%%)", length(het.dia.patient)/length(levels(m.relapsing.dia$patient))*100), sprintf("homogenious (%.0f%%)", 100-length(het.dia.patient)/length(levels(m.relapsing.dia$patient))*100))
m.relapsing.dia$group <- ifelse(m.relapsing.dia$patient %in% het.dia.patient, "heterogeneous", "homogeneous")
m.relapsing.dia$label <- ifelse(paste0(m.relapsing.dia$patient, ":", m.relapsing.dia$mut) %in% m.relapsing.cons, "symbol(\"\\267\")", NA)

plot.dia <- ggplot(data=m.relapsing.dia, aes(x=patient.label, y=frequency.norm, fill=mut, order=-frequency.norm)) + 
		facet_grid(.~group, scale="free_x", space = "free_x") +
		geom_bar(stat="identity", width=0.9, colour="white") +
		geom_text(aes(label=label, y=midpoint), size=8, vjust=1, colour="white", parse=TRUE) +
		scale_fill_manual(values = cols, breaks=c("KRAS:12:25398284:C>T", "NRAS:1:115258744:C>T", "PTPN11:12:112888210:G>A", "FLT3:13:28599081:C>A", "JAK1:1:65305426:G>C", "JAK2:9:5073749:G>A", "JAK3:19:17945970:G>A", "CRLF2:X:1314966:A>C"), labels=c("KRAS", "NRAS", "PTPN11", "FLT3", "JAK1", "JAK2", "JAK3", "CRLF2")) + 
		scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks=c(0.25, 0.5, 0.75, 1)) +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.text.y = element_text(size=13), axis.title.y = element_text(size=18, vjust=0.1), axis.title.x = element_blank(), 
				legend.title = element_blank(), legend.text=element_text(size=13), legend.key.size=unit(0.7, "cm"), legend.justification=c(1,1), legend.position=c(1,1), legend.background=element_rect(colour="black"),
				strip.text.x = element_text(size=15), 
				panel.grid.major.x = element_blank(),
				panel.grid.minor.y = element_blank(),
				plot.title=element_text(size=20, face="bold", vjust=1), plot.margin=unit(c(0.5,1,0,1), "cm")) +
#		geom_text(aes(label = mut.short), size = 3, hjust = 0.5, vjust = 1.2, position = "stack", colour="white") +
		ylab("Adjusted AF diagnosis")	+ 
		ggtitle("Ras/Jak signaling mutations at diagnosis")

#---
# RELAPSING, REL
#---

# reorder patients by cumulative frequency
sorted <- aggregate(frequency.norm~patient, data=m.relapsing.rel, FUN=sum)
sorted <- sorted[order(sorted$frequency.norm, decreasing=T), c("patient", "frequency.norm")]
m.relapsing.rel$patient <- factor(m.relapsing.rel$patient, sorted$patient)

# find heterogenious samples
het.rel <- aggregate(frequency.norm~patient, data=m.relapsing.rel, FUN=length)
het.rel.patient <- het.rel[het.rel$frequency.norm>1,"patient"]
#m.relapsing.rel$group <- ifelse(m.relapsing.rel$patient %in% het.rel.patient, sprintf("heterogenious (%.0f%%)", length(het.rel.patient)/length(levels(m.relapsing.rel$patient))*100), sprintf("homogenious (%.0f%%)", 100-length(het.rel.patient)/length(levels(m.relapsing.rel$patient))*100))
m.relapsing.rel$group <- ifelse(m.relapsing.rel$patient %in% het.rel.patient, "heterogeneous", "homogeneous")
m.relapsing.rel$label <- ifelse(paste0(m.relapsing.rel$patient, ":", m.relapsing.rel$mut) %in% m.relapsing.cons, "symbol(\"\\267\")", NA)

plot.rel <- ggplot(data=m.relapsing.rel, aes(x=patient, y=frequency.norm, fill=mut, order=-frequency.norm)) + 
				facet_grid(.~group, scale="free_x", space = "free_x") +
				geom_bar(stat="identity", width=0.9, colour="white") +
				geom_text(aes(label=label, y=midpoint), size=8, vjust=1, colour="white", parse=TRUE) +
				scale_fill_manual(values = cols, name="Mutation", guide=FALSE) + 
				scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks=c(0.25, 0.5, 0.75, 1)) +
				theme_bw() + 
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.text.y = element_text(size=13), axis.title.x = element_text(size=20, vjust=0.2), axis.title.y = element_text(size=18, vjust=0.1), 
						legend.key.size=unit(0.35, "cm"),
						strip.text.x = element_text(size=15), 
						panel.grid.major.x = element_blank(),
						panel.grid.minor.y = element_blank(),
						plot.title=element_text(size=20, face="bold", vjust=1), plot.margin=unit(c(1,1,0.3,1), "cm")) +
#				geom_text(aes(label = mut.short), size = 3, hjust = 0.5, vjust = 1.2, position = "stack", colour="white") +
				ylab("Adjusted AF relapse") + 
				ggtitle("Ras/Jak signaling mutations at relapse")
		
grid.arrange(plot.dia, plot.rel, nrow=2)
dev.off()

#---
# RELAPSING, MATCHED, DIA+REL
#---

# reorder patients by cumulative frequency
m.relapsing.combined <- m.relapsing.dia[m.relapsing.dia$patient %in% cases.relapsing.matched,]
m.relapsing.combined <- rbind(m.relapsing.combined, m.relapsing.rel[m.relapsing.rel$patient %in% cases.relapsing.matched,])

# reorder patients by group and cumulative frequency
sorted <- aggregate(frequency.norm~patient+sample+group, data=m.relapsing.combined, FUN=sum)
sorted$sample <- as.character(sorted$sample)
sorted$sample[sorted$sample=="relapse"] <- "0relapse"
sorted$group[sorted$group=="homogeneous"] <- "0homogeneous"
sorted <- sorted[order(sorted$sample, sorted$group, sorted$frequency.norm, decreasing=T),]
m.relapsing.combined$patient <- factor(m.relapsing.combined$patient, unique(as.character(sorted$patient)))

pdf("/mnt/projects/p2ry8-crlf2/results/figures/signaling-barplot.matched.pdf", width=14, height=9)
plot.dia <- ggplot(data=m.relapsing.combined, aes(x=patient, y=frequency.norm, fill=mut, order=-frequency.norm)) + 
		facet_grid(sample~.) +
		geom_bar(stat="identity", width=0.9, colour="white") +
		geom_text(aes(label=label, y=midpoint), size=8, vjust=1, colour="white", parse=TRUE) +
		scale_fill_manual(values = cols, breaks=c("KRAS:12:25398284:C>T", "NRAS:1:115258744:C>T", "PTPN11:12:112888210:G>A", "FLT3:13:28599081:C>A", "JAK1:1:65305426:G>C", "JAK2:9:5078360:A>G", "JAK3:19:17945970:G>A", "CRLF2:X:1314966:A>C"), labels=c("KRAS", "NRAS", "PTPN11", "FLT3", "JAK1", "JAK2", "JAK3", "CRLF2")) + 
		scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks=c(0.25, 0.5, 0.75, 1)) +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.text.y = element_text(size=13), axis.title.y = element_text(size=18, vjust=0.1), axis.title.x = element_blank(), 
				legend.title = element_blank(), legend.text=element_text(size=13), legend.key.size=unit(0.7, "cm"), legend.justification=c(1,1), legend.position=c(1,1), legend.background=element_rect(colour="black"),
				strip.text.x = element_text(size=15), 
				strip.text.y = element_text(size=15), 
				panel.grid.major.x = element_blank(),
				panel.grid.minor.y = element_blank(),
				plot.title=element_text(size=20, face="bold", vjust=1), plot.margin=unit(c(0.5,1,0,1), "cm")) +
#		geom_text(aes(label = mut.short), size = 3, hjust = 0.5, vjust = 1.2, position = "stack", colour="white") +
		ylab("Adjusted AF") + 
		ggtitle("Ras/Jak signaling mutations at diagnosis and relapse (matched cases)")
print(plot.dia)
dev.off()

#---
# NON-RELAPSING, DIA
#---

# reorder patients by cumulative frequency
sorted <- aggregate(frequency.norm~patient, data=m.nonrelapsing, FUN=sum)
sorted <- sorted[order(sorted$frequency.norm, decreasing=T), c("patient", "frequency.norm")]
m.nonrelapsing$patient <- factor(m.nonrelapsing$patient, sorted$patient)

# find heterogenious samples
het <- aggregate(frequency.norm~patient, data=m.nonrelapsing, FUN=length)
het <- het[het$frequency.norm>1,"patient"]
m.nonrelapsing$group <- ifelse(m.nonrelapsing$patient %in% het, "heterogeneous", "homogeneous")

pdf("/mnt/projects/p2ry8-crlf2/results/figures/signaling-barplot.non-relapsing.pdf", width=10, height=6)
print(ggplot(data=m.nonrelapsing, aes(x=patient, y=frequency.norm, fill=mut, order=-frequency.norm)) +
				facet_grid(.~group, scale="free_x", space = "free_x") +
				geom_bar(stat="identity", width=0.9, colour="white") +
				scale_fill_manual(values = cols, breaks=c("KRAS:12:25398284:C>T", "NRAS:1:115258747:C>T", "PTPN11:12:112888210:G>A", "JAK1:1:65305426:G>C", "JAK2:9:5078395:C>G", "JAK3:19:17945969:C>T", "CRLF2:X:1314966:A>C"), labels=c("KRAS", "NRAS", "PTPN11", "JAK1", "JAK2", "JAK3", "CRLF2")) + 
				#scale_fill_manual(values = cols, name="Mutation") + 
				scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks=c(0.25, 0.5, 0.75, 1)) +
				theme_bw() + 
				theme(axis.text.x = element_text(angle = 90, hjust = 1),
					  legend.key.size=unit(0.35, "cm"),
					  panel.grid.major.x = element_blank(),
					  panel.grid.minor.y = element_blank()) +
#				geom_text(aes(label = mut.short), size = 3, hjust = 0.5, vjust = 1.2, position = "stack", colour="white") +
				ylab("Adjusted AF") +
				ggtitle("Ras/Jak signaling mutations at diagnosis (non-relapsing cohort)"))
dev.off()

# compute significance for difference in homogeneous/heterogenious samples b/w diagnosis/relapse
cont <- matrix(c(sum(het.dia$frequency.norm>1), sum(het.dia$frequency.norm==1), sum(het.rel$frequency.norm>1), sum(het.rel$frequency.norm==1)), nrow=2, byrow=T)
print(cont)
print(sprintf("P-value homogeneous/heterogenious b/w diagnosis/relapse (Fisher's exact test): %g", fisher.test(cont)$p.value))

# compute significance for difference in ras/jak mutations b/w diagnosis/relapse
t <- data.frame(patient=rep(cases.relapsing.matched,2), sample=c(rep("diagnosis", length(cases.relapsing.matched)), rep("relapse", length(cases.relapsing.matched))))
t <- merge(t, aggregate(gene~patient+sample, m, paste), all.x=T)
ras <- grepl("(KRAS|NRAS|PTPN11|FLT3)", t$gene)
jak <- grepl("(JAK1|JAK2|JAK3|CRLF2)", t$gene)
t$sigstat <- NA
t$sigstat[is.na(t$gene)] <- "none"
t$sigstat[ras & !jak] <- "RAS"
t$sigstat[!ras & jak] <- "JAK"
t$sigstat[ras & jak] <- "both"
cont <- table(t$sample, t$sigstat)
print(cont)
print(sprintf("P-value ras/jak distribution b/w diagnosis/relapse (Fisher's exact test): %g", fisher.test(cont)$p.value))
