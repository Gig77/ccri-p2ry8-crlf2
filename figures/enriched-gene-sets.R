keep <- c("NAME", "CATEGORY", "LINKOUT", "SIZE", "NES", "NOM.p.val", "FDR.q.val")
merge.by <- c("NAME", "CATEGORY", "LINKOUT", "SIZE")
q.cutoff <- 0.1

up <- read.delim("/mnt/projects/ikaros/results/anduril/execute/gseaMSigDB_IKNp_vs_IKCp/enrichedUp.csv", stringsAsFactors = F)[,keep]
names(up)[5:7] <- c("NES.IKNpVsIKCp", "p.IKNpVsIKCp", "q.IKNpVsIKCp")
dn <- read.delim("/mnt/projects/ikaros/results/anduril/execute/gseaMSigDB_IKNp_vs_IKCp/enrichedDown.csv", stringsAsFactors = F)[,keep]
names(dn)[5:7] <- c("NES.IKNpVsIKCp", "p.IKNpVsIKCp", "q.IKNpVsIKCp")
gs <- merge(up, dn, by=c(merge.by, "NES.IKNpVsIKCp", "p.IKNpVsIKCp", "q.IKNpVsIKCp"), all=T)

up <- read.delim("/mnt/projects/ikaros/results/anduril/execute/gseaMSigDB_IKMp_vs_IKCp/enrichedUp.csv", stringsAsFactors = F)[,keep]
names(up)[5:7] <- c("NES.IKMpVsIKCp", "p.IKMpVsIKCp", "q.IKMpVsIKCp")
dn <- read.delim("/mnt/projects/ikaros/results/anduril/execute/gseaMSigDB_IKMp_vs_IKCp/enrichedDown.csv", stringsAsFactors = F)[,keep]
names(dn)[5:7] <- c("NES.IKMpVsIKCp", "p.IKMpVsIKCp", "q.IKMpVsIKCp")
gs <- merge(gs, merge(up, dn, by=c(merge.by, "NES.IKMpVsIKCp", "p.IKMpVsIKCp", "q.IKMpVsIKCp"), all=T), by=merge.by, all=T)

up <- read.delim("/mnt/projects/ikaros/results/anduril/execute/gseaMSigDB_IKD_vs_IKCp/enrichedUp.csv", stringsAsFactors = F)[,keep]
names(up)[5:7] <- c("NES.IKDVsIKCp", "p.IKDVsIKCp", "q.IKDVsIKCp")
dn <- read.delim("/mnt/projects/ikaros/results/anduril/execute/gseaMSigDB_IKD_vs_IKCp/enrichedDown.csv", stringsAsFactors = F)[,keep]
names(dn)[5:7] <- c("NES.IKDVsIKCp", "p.IKDVsIKCp", "q.IKDVsIKCp")
gs <- merge(gs, merge(up, dn, by=c(merge.by, "NES.IKDVsIKCp", "p.IKDVsIKCp", "q.IKDVsIKCp"), all=T), by=merge.by, all=T)

gs <- gs[apply(gs[,grep("q.", names(gs))], 1, function(x) min(x, na.rm=T) <= q.cutoff) > 0,]
write.table(gs, "/mnt/projects/p2ry8-crlf2/results/figures/enriched-gene-sets.tsv", row.names = F, col.names = T, sep = "\t", quote = F, na = "")
