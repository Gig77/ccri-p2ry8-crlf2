t <- read.delim("/mnt/projects/p2ry8-crlf2/results/exomeCopy/allpatients.filtered-segments.exomeCopy.tsv", stringsAsFactors=F, check.names = F)

names(t)[names(t)=="seqnames"] <- "chr"
names(t)[names(t)=="sample.name"] <- "sample"
names(t)[names(t)=="nranges"] <- "nexons"
t <- t[,c("sample", "chr", "start", "end", "width", "nexons", "targeted.bp", "copy.count", "log.odds", "genes")]

# excluded samples
t <- t[!t$sample %in% c("LU3D", "SN18D", "564D", "460D", "545D", "MA5D", "957D", "KE17247D", "BJ17183D"),]
t <- t[!t$sample %in% c("LU3R", "SN18R", "564R", "460R", "545R", "MA5R", "957R"),]
t <- t[!t$sample %in% c("737R2"),]

# exclude xenografts
t <- t[!grepl("^m\\d", t$sample),]

# exclude low-quality samples
t <- t[!t$sample %in% c("715D", "715R", "GI8R"),]

# rename some samples
t$sample[t$sample=="108R2"] <- "108RR"
t$sample[t$sample=="715R3"] <- "715RR"
t$sample[t$sample=="737R3"] <- "737RR"
t$sample[t$sample=="AL9890R2"] <- "AL9890RR"
t$sample[t$sample=="S23R3"] <- "S23RR"

# exclude remission samples
t <- t[!grepl("C$", t$sample),]

# sort
t <- t[order(t$sample, t$chr, t$start),]

# write
write.table(t, file="/mnt/projects/p2ry8-crlf2/results/exomeCopy/filtered-cnvs.paper.tsv.part", sep="\t", col.names=T, row.names=F, quote=F, na="")
