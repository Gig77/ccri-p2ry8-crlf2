# primary samples
var.prim <- read.delim("/mnt/projects/p2ry8-crlf2/results/filtered-variants.cosmic.tsv", check.names=F, stringsAsFactor=F)
var.prim$sample <- gsub("rem_", "", var.prim$sample)

# xenografts
var.xeno <- read.delim("/mnt/projects/p2ry8-crlf2/results/filtered-variants.xenografts.cosmic.tsv", check.names=F, stringsAsFactor=F)
var.xeno$sample[var.xeno$patient %in% c("m1977-G-dia", "m1966-Y-dia", "m1035-108-dia", "m252-379-dia", "m1041-737-dia", "m247-833-dia", "m1059-92-dia", "m1037-839-dia", "m248-841-dia")] <- "dia-xeno"
var.xeno$sample[var.xeno$patient %in% c("m1963-545-rel", "m1957-715-rel", "m1967-Y-rel", "m1060-108-rel", "m1069-737-rel")] <- "rel-xeno"
var.xeno$sample[var.xeno$patient %in% c("m1964-545-rel")] <- "rel-xeno2"
var.xeno$patient <- gsub("(.*)-(.*)-(.*)", "\\2", var.xeno$patient)

# combine and filter
var <- rbind(var.prim, var.xeno)
var <- var[var$status!="REJECT",]

# merge
fields.common <- c("patient", "cohort", "var_type", "chr", "pos", "ref", "alt", "gene", "effect", "non_silent", "deleterious", "aa_change", "Polyphen2", "SIFT", "SiPhy", "cosmic_hits_nt", "cosmic_hits_aa")
fields.sample <- c("dp_leu_tot", "dp_leu_ref", "dp_leu_var", "freq_leu")

m <- var[var$sample=="dia", c(fields.common, fields.sample)]
names(m)[names(m) %in% fields.sample] <- c("dp_tot.dia", "dp_ref.dia", "dp_var.dia", "freq.dia")

m <- merge(m, var[var$sample=="rel",c(fields.common, fields.sample)], by=fields.common, all=TRUE)
names(m)[names(m) %in% fields.sample] <- c("dp_tot.rel", "dp_ref.rel", "dp_var.rel", "freq.rel")

m <- merge(m, var[var$sample=="rel2",c(fields.common, fields.sample)], by=fields.common, all=TRUE)
names(m)[names(m) %in% fields.sample] <- c("dp_tot.rel2", "dp_ref.rel2", "dp_var.rel2", "freq.rel2")

m <- merge(m, var[var$sample=="rel3",c(fields.common, fields.sample)], by=fields.common, all=TRUE)
names(m)[names(m) %in% fields.sample] <- c("dp_tot.rel3", "dp_ref.rel3", "dp_var.rel3", "freq.rel3")

m <- merge(m, var[var$sample=="dia-xeno",c(fields.common, fields.sample)], by=fields.common, all.x=TRUE)
names(m)[names(m) %in% fields.sample] <- c("dp_tot.dia-xeno", "dp_ref.dia-xeno", "dp_var.dia-xeno", "freq.dia-xeno")

m <- merge(m, var[var$sample=="rel-xeno",c(fields.common, fields.sample)], by=fields.common, all.x=TRUE)
names(m)[names(m) %in% fields.sample] <- c("dp_tot.rel-xeno", "dp_ref.rel-xeno", "dp_var.rel-xeno", "freq.rel-xeno")

m <- merge(m, var[var$sample=="rel-xeno2",c(fields.common, fields.sample)], by=fields.common, all.x=TRUE)
names(m)[names(m) %in% fields.sample] <- c("dp_tot.rel-xeno2", "dp_ref.rel-xeno2", "dp_var.rel-xeno2", "freq.rel-xeno2")

write.table(m, file="/mnt/projects/p2ry8-crlf2/results/filtered-variants.cosmic.merged.tsv", row.names=F, sep="\t", quote=F, na="")
