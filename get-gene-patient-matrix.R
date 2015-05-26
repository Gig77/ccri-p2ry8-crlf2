library(reshape)
warnings()

t <- read.delim("/mnt/projects/p2ry8-crlf2/results/filtered-variants.tsv", stringsAsFactors=F)
t <- t[t$status != "REJECT" & t$freq_leu >= 0.1 & t$non_silent==1,]

t$sample[t$sample=="rem_dia"] <- "dia"
t$sample[t$sample=="rem_rel"] <- "rel"
t$sample[t$sample=="rem_rel2"] <- "rel2"
t$sample[t$sample=="rem_rel3"] <- "rel3"

# determine relapsing patients for which we have both diagnosis and relapse samples
patients.relapsing <- unique(t$patient[t$sample=="dia"])[unique(t$patient[t$sample=="dia"]) %in% unique(t$patient[t$sample=="rel"])]

# flag conserved variants
conserved <- merge(t[t$sample=="dia", c("patient", "chr", "pos", "ref", "alt")], t[t$sample=="rel", c("patient", "chr", "pos", "ref", "alt")])
conserved$conserved = "yes"
t <- merge(t, conserved, all.x=T)
t$conserved[is.na(t$conserved)] <- ifelse(t$patient[is.na(t$conserved)] %in% patients.relapsing, "no", NA)

# format cell value
t$cv <- sprintf("%d%s%s", as.integer(t$freq_leu*100), ifelse(!is.na(t$deleterious) & t$deleterious == "yes", "d", ""), ifelse(!is.na(t$conserved) & t$conserved == "yes", "c", ""))

# reformat table
m <- cast(t, formula = gene ~ patient + sample, value="cv", fun.aggregate=function(x) { paste(x, collapse=";") }, fill="")

# order columns
m <- m[,c(1,order(grepl("dia", names(m)[-1]), m[m$gene=="JAK2",-1] != "", m[m$gene=="JAK3",-1] != "", m[m$gene=="KRAS",-1] != "", m[m$gene=="NRAS",-1] != "", decreasing=T)+1)]

# insert sums
num_cases_rel <- apply(m[,-1], 1, function(x) { sum(x != "" & grepl("rel", names(m)[-1])) })
num_cases_dia <- apply(m[,-1], 1, function(x) { sum(x != "" & grepl("dia", names(m)[-1])) })
m <- data.frame(m["gene"], num_cases_dia, m[grep("dia", names(m))], num_cases_rel, m[grep("rel", names(m))], check.names=F)

# order rows
topgenes <- c("JAK2", "JAK3", "JAK1", "KRAS", "NRAS", "PTPN11", "FLT3", "SETD2", "CREBBP", "ETV6", "PAX5", "IKZF1", "EP300", "BRCA2", "MET", "MSH6", "CDKN2A", "IL7R")
m <- m[unlist(c(sapply(topgenes, function(x) { which(x==m$gene) }), which(!m$gene %in% topgenes))),]

# indicate diagnosis samples from relapsing cohort with "_r" suffix in sample name
names(m)[names(m) %in% paste0(patients.relapsing, "_dia")] <- paste0(names(m)[names(m) %in% paste0(patients.relapsing, "_dia")], "_r")

# write output
write.table(m, file="/mnt/projects/p2ry8-crlf2/results/gene-patient-matrix.tsv.part", row.names=F, col.names=T, sep="\t")
