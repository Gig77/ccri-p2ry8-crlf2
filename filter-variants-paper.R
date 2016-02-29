t <- read.delim("/mnt/projects/p2ry8-crlf2/results/filtered-variants.cosmic.tsv", stringsAsFactors=F, check.names = F)

t <- t[t$status != "REJECT" & t$freq_leu >= 0.1 & !grepl("^GL", t$chr),]
t$status[t$status %in% c(".", "MISSED")] <- "PASS"

names(t)[names(t)=="sample"] <- "timepoint"

# excluded samples
t <- t[!t$patient %in% c("LU3", "SN18", "564", "460", "545", "MA5", "957"),]
t <- t[t$patient != "KE17247" | t$timepoint != "rem_rel",]
t <- t[t$patient != "BJ17183" | t$timepoint != "rem_rel",]
t <- t[t$patient != "737" | t$timepoint != "rem_rel2",]  # only 8% blasts

# rename timepoints
t$timepoint[t$timepoint=="rem_dia"] <- "D"
t$timepoint[t$timepoint=="rem_rel"] <- "R"
t$timepoint[t$timepoint %in% c("rem_rel2", "rem_rel3")] <- "RR"

# exclude some columns
t <- t[,!names(t) %in% c("rejected_because", "rem_samples", "cosmic_hits_leu_nt", "cosmic_hits_leu_aa")]

# sort
t <- t[order(t$patient, t$timepoint, t$chr, t$pos),]

# write
write.table(t, file="/mnt/projects/p2ry8-crlf2/results/filtered-variants.paper.tsv.part", sep="\t", col.names=T, row.names=F, quote=F, na="")
