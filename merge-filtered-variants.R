var <- read.delim("/mnt/projects/p2ry8-crlf2/results/filtered-variants.cosmic.tsv", check.names=F, stringsAsFactor=F)
var.pass <- var[var$status!="REJECT",]
dia <- var.pass[var.pass$sample=="rem_dia",]
rel1 <- var.pass[var.pass$sample=="rem_rel",]
rel2 <- var.pass[var.pass$sample=="rem_rel2",]
m <- merge(dia, rel1, by=c("patient", "cohort", "status", "var_type", "chr", "pos", "dbSNP", "ref", "alt", "gene", "add_genes", "aa_change", "impact", "effect", "snpeff_effect", "non_silent", "deleterious", "exons", "dp_rem_tot", "dp_rem_ref", "dp_rem_var", "freq_rem", "cosmic_hits_nt", "cosmic_hits_aa", "cosmic_hits_leu_nt", "cosmic_hits_leu_aa"), all=T, suffixes=c(".dia", ".rel1"))
m2 <- merge(m, rel2, by=c("patient", "cohort", "status", "var_type", "chr", "pos", "dbSNP", "ref", "alt", "gene", "add_genes", "aa_change", "impact", "effect", "snpeff_effect", "non_silent", "deleterious", "exons", "dp_rem_tot", "dp_rem_ref", "dp_rem_var", "freq_rem", "cosmic_hits_nt", "cosmic_hits_aa", "cosmic_hits_leu_nt", "cosmic_hits_leu_aa"), all=T, suffixes=c("", ".rel2"))
write.table(m2, file="filtered-variants.cosmic.merged.tsv", row.names=F, sep="\t", quote=F, na="")

