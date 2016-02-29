library("RColorBrewer")
warnings()

# TABLE: filtered-variants.tsv
# read and filter input data
t <- read.delim("/mnt/projects/p2ry8-crlf2/results/filtered-variants.tsv", stringsAsFactors=F)
t <- t[t$status != "REJECT" & t$freq_leu >= 0.1 & t$non_silent==1,]

# excluded samples
t <- t[!t$patient %in% c("LU3", "SN18", "564", "460", "545", "MA5", "957"),]
t <- t[t$patient != "KE17247" | t$sample != "rem_rel",]
t <- t[t$patient != "BJ17183" | t$sample != "rem_rel",]

# join diagnosis and relapse variants
t.dia <- t[t$sample == "rem_dia", c("patient", "cohort", "chr", "pos", "ref", "alt", "freq_leu")]
names(t.dia)[7] <- "dia"
t.rel <- t[t$sample == "rem_rel", c("patient", "cohort", "chr", "pos", "ref", "alt", "freq_leu")]
names(t.rel)[7] <- "rel"
m <- merge(t.dia, t.rel, all.x=T, all.y=T)

patients <- unique(m$patient)
patients.relapsing <- patients[patients %in% t$patient[t$sample=="rem_rel"]]

bardata <- matrix(rep(0, 3*length(patients)), nrow=length(patients))
colnames(bardata) <- c("relapse-specific", "diagnosis-specific", "conserved")
rownames(bardata) <- patients

# bar data
for(p in patients) {
	m.p <- m[m$patient==p,]
	bardata[p,1] <- ifelse(p %in% patients.relapsing, sum(is.na(m.p$dia) & !is.na(m.p$rel)), NA)
	bardata[p,2] <- sum(!is.na(m.p$dia) & is.na(m.p$rel))
	bardata[p,3] <- ifelse(p %in% patients.relapsing, sum(!is.na(m.p$dia) & !is.na(m.p$rel)), NA)
}
bardata <- bardata[order(rowSums(bardata)),]
bardata.relapsing <- bardata[!is.na(bardata[,1]),]

# boxplot data
mut <- as.data.frame(bardata)
mut$dia <- ifelse(is.na(mut$conserved), mut$'diagnosis-specific', mut$'diagnosis-specific' + mut$conserved)
mut$rel <- ifelse(is.na(mut$conserved), NA, mut$'relapse-specific' + mut$conserved)
#test <- kruskal.test(list(mut$dia, mut$rel))

# plot
make_plot <- function() {
	par(mfrow=c(2,2), mar=c(7,5,1,1))
	
	# boxplot (include outliers)
	boxplot(mut$dia, mut$rel, log="y", xlab="", ylab="# non-silent mutations", na.action=na.exclude, outline=F, cex.axis=1.4, xaxt="n", cex.lab=1.5, ylim=c(5,1000))
	axis(1, at=c(1,2), cex.axis=1.5, labels=c(sprintf("diagnosis (n=%d)", sum(!is.na(mut$dia))), sprintf("relapse (n=%d)", sum(!is.na(mut$rel)))))
	stripchart(list(mut$dia, mut$rel), method="jitter", vertical=T, pch=19, col=c("red", "blue"), add=T)
	#print(sprintf("P-value Kruskal-Wallis: %.2g", test$p.value))

	# boxplot (exclude outliers)
	boxplot(mut$dia, mut$rel, xlab="", ylab="# non-silent mutations", na.action=na.exclude, outline=F, cex.axis=1.5, xaxt="n", cex.lab=1.5)
	axis(1, at=c(1,2), cex.axis=1.5, labels=c(sprintf("diagnosis (n=%d)", sum(!is.na(mut$dia))), sprintf("relapse (n=%d)", sum(!is.na(mut$rel)))))
	stripchart(list(mut$dia, mut$rel), method="jitter", vertical=T, pch=19, col=c("red", "blue"), add=T)
	
	# barplot (include outliers)
	barplot(t(bardata.relapsing), beside=FALSE, col=brewer.pal(6, "Set2"), las=2, yaxt="n", ylim=c(0,max(rowSums(bardata.relapsing))), cex.axis=1.5, cex.names=1.5, cex.lab=1.5)
	abline(h=seq(0, max(rowSums(bardata.relapsing)), by=100), col="gray")
	box()
	barplot(t(bardata.relapsing), beside=FALSE, col=brewer.pal(6, "Set2"), ylab="# non-silent mutations", xaxt="n", ylim=c(0,max(rowSums(bardata.relapsing))), cex.axis=1.5, cex.names=1.5, cex.lab=1.5, add=T)
	legend(x=1, y=800, bg="white", rev(colnames(bardata.relapsing)), fill=rev(brewer.pal(3, "Set2")), cex=1.5)

	# barplot (cut-off outliers)
	maxy <- 200
	barplot(t(bardata.relapsing), beside=FALSE, col=brewer.pal(6, "Set2"), las=2, yaxt="n", ylim=c(0,maxy), cex.axis=1.5, cex.names=1.5, cex.lab=1.5)
	abline(h=seq(0, maxy, by=10), col="gray")
	box()
	barplot(t(bardata.relapsing), beside=FALSE, col=brewer.pal(6, "Set2"), ylab="# non-silent mutations", xaxt="n", ylim=c(0,maxy), cex.axis=1.5, cex.names=1.5, cex.lab=1.5, add=T)
	legend(x=1, y=190, bg="white", rev(colnames(bardata.relapsing)), fill=rev(brewer.pal(3, "Set2")), cex=1.5)
}

# plot PNG
png("/mnt/projects/p2ry8-crlf2/results/figures/mutations-per-patient.png", width=2800, height=2800, pointsize=40)
make_plot()
dev.off()

# plot PDF
pdf("/mnt/projects/p2ry8-crlf2/results/figures/mutations-per-patient.pdf", width=13, height=13)
make_plot()
dev.off()
