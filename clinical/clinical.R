options(warn=1)
library(gridExtra)
library(vcd)

rm(list=ls())
source("~/p2ry8-crlf2/scripts/clinical/test_pairwise.R")

set.seed(22)

patients.relapsing.matched <- c("108", "92", "715", "737", "839", "B36", "BB16", "DL2", "GI8", "GI13", "HV57", "HV80", "N7", "S23", "DS10898", "VS14645", "SE15285", "KE17247", "AL9890", "GL11356")
patients.relapsing.diaonly <- c("1060", "BJ17183")
patients.relapsing <- c(patients.relapsing.matched, patients.relapsing.diaonly)
patients.nonrelapsing.withmutdata <- c("242", "360", "365", "379", "400", "506", "769", "802", "833", "887", "841", "903", "948", "961", "1066", "1089", "HW11537", "KT14158", "TL14516")
patients.nonrelapsing.nomutdata <- c("5755", "7839", "11898", "14197", "6603", "7118", "11536", "7361", "13906", "4558", "4868")
patients.nonrelapsing <- c(patients.nonrelapsing.withmutdata, patients.nonrelapsing.nomutdata)
patients <- c(patients.relapsing, patients.nonrelapsing)
patients.excl <- c("MA5", "BJ14367", "LU3", "SN18", "460", "545", "564", "957")

c <- read.delim("~/p2ry8-crlf2/results/clinical/Clinical data_P-C_v8.txt", strip.white=TRUE, na.strings=c("", "NA", "-", " -", " N.A.", "N.A", "N.A.", "N.A.?", "n/a", "n/d", "n.d.", " ", "early (CNS)"))

#----
# CLEAN UP DATA
#----

c$patient_id <- as.character(c$patient_id)
c$patient_id[c$patient_id=="N 7"] <- "N7"
c$patient_id <- as.factor(c$patient_id)
sprintf("Patient not in clinical data sheet: %s", paste(patients[!patients %in% c$patient_id], sep=","))
sprintf("The following patient is in the clinical data but without mutation data: %s", paste(c$patient_id[!c$patient_id %in% c(patients, patients.excl)], sep=","))
sprintf("Patient deliberately excluded from analysis: %s", paste(patients.excl, sep=","))
c <- c[c$patient_id %in% patients,]
c$cohort[c$cohort=="relapsing EM"] <- "relapsing"
c$cohort <- as.factor(as.character(c$cohort))
c$blasts_dia[c$blasts_dia=="73 pB"] <- 73
c$blasts_dia <- as.numeric(as.character(c$blasts_dia))
c$mrd_risk_dia[c$mrd_risk_dia=="SER"] <- "IR"
c$mrd_risk_dia[c$mrd_risk_dia=="LR"] <- "SR"
c$mrd_risk_dia <- as.factor(as.character(c$mrd_risk_dia))
c$rel_timepoint[c$rel_timepoint=="very early"] <- "early"
c$rel_timepoint <- as.factor(as.character(c$rel_timepoint))
c$outcome_after.1st.rel[c$patient_id=="1060"] <- NA
c$outcome_after.1st.rel <- as.factor(as.character(c$outcome_after.1st.rel))
c$endpoint[c$patient_id=="1060"] <- NA
c$second_rem_months[c$patient_id=="1060"] <- NA
names(c)[names(c)=="blasts_rel..1st."] <- "blasts_rel"
c$blasts_rel[c$blasts_rel=="isoliertes ZNS Rezidiv"] <- NA
c$blasts_rel <- as.numeric(as.character(c$blasts_rel))
c$rel_site[c$rel_site=="BM combined (testes)"] <- "BM combined"
c$rel_site <- as.factor(as.character(c$rel_site))
c$mrd_level_rel[c$mrd_level_rel=="poor (78%Bl.)"] <- "poor"
c$mrd_level_rel <- as.factor(as.character(c$mrd_level_rel))
c$rel_protocol <- as.factor(toupper(as.character(c$rel_protocol)))
c$iAMP21 <- as.logical(ifelse(is.na(c$iAMP21), 0, 1))
c$DS <- as.logical(ifelse(is.na(c$DS), 0, 1))
c$HD <- as.logical(ifelse(is.na(c$HD), 0, 1))

#----
# READ MUTATION DATA
#----

m <- read.delim("~/p2ry8-crlf2/results/filtered-variants.cosmic.tsv")
m <- m[m$status!="REJECT" & m$non_silent==T & m$freq_leu >= 0.1,]
m$patient <- as.character(m$patient)
m.merged <- merge(m[m$sample=="rem_dia", c("patient", "chr", "pos", "ref", "alt", "freq_leu", "gene")], m[m$sample=="rem_rel",c("patient", "chr", "pos", "ref", "alt", "freq_leu", "gene")], by=c("patient", "chr", "pos", "ref", "alt", "gene"), all=T, suffixes=c(".dia", ".rel")) 

#cn <- read.delim("~/p2ry8-crlf2/results/exomeCopy/allpatients.filtered-segments.exomeCopy.tsv")
#cn$patient <- unlist(strsplit(as.character(cn$sample.name), "(C|D|R\\d?)$", perl=T))
#cn$sample <- NA
#cn$sample[grepl("D$", cn$sample.name, perl=T)] <- "dia"
#cn$sample[grepl("C$", cn$sample.name, perl=T)] <- "rem"
#cn$sample[grepl("R\\d?$", cn$sample.name, perl=T)] <- "rel"
#cn <- cn[cn$sample != "rem",]
#cn <- cn[cn$patient %in% c("564", "715")] # exclude noisy samples

ikzf1.del.dia <- c("108", "DL2", "737", "841", "92", "BJ17183", "DS10898", "HV57", "HV80", "KT14158") 
ikzf1.del.rel <- c("108", "DL2", "737", "92", "B36", "DS10898", "GI8", "GI13", "HV57", "HV80", "N7")

ikzf2.del.dia <- c("")
ikzf2.del.rel <- c("HV57", "KE17247")

ikzf3.del.dia <- c("1089")
ikzf3.del.rel <- c("")

cdkn2a.del.dia <- c("1060", "360", "379", "400", "506", "564", "737", "769", "887", "957", "B36", "BB16", "LU3", "S23", "SN18", "14197", "6603", "7118", "11536", "7361", "13906")
cdkn2a.del.rel <- c("108", "545", "564", "737", "B36", "BB16", "GI8", "HV80", "SE15285", "SN18")

pax5.del.dia <- c("1060", "365", "379", "400", "545", "564", "737", "833", "948", "957", "BJ17183", "HV80", "KT14158", "5755", "7839", "6603", "7361", "903") 
pax5.del.rel <- c("545", "564", "737", "DL2", "HV80", "N7", "VS14645")

ebf1.del.dia <- c("839", "GL11356")
ebf1.del.rel <- c("839", "GL11356")

etv6.del.dia <- c("N7", "242", "957", "11898", "14197")
etv6.del.rel <- c("N7", "715", "GI13")

setd2.del.dia <- c("903", "VS14645", "KE17247", "AL9890")
setd2.del.rel <- c("VS14645", "DS10898")

p2ry8.sub <- c("545", "460", "564", "SN18", "LU3", "957")
p2ry8.cons <- c("DL2", "108", "GI8", "HV80", "N7", "DS10898", "VS14645", "BB16", "SE15285", "715", "839", "AL9890", "GL11356")

chr21.gain <- c("108", "DL2", "N7", "DS10898", "SE15285", "GI13", "HV57", "737", "841", "365", "903", "400", "HW11537", "TL14516", "887", "769", "506", "1089", "802", "715", "545", "460", "1066", "961", "GI8", "1060", "VS14645", "360", "564", "5755", "7839", "13906", "4558", "957", "AL9890", "GL11356")

sex.chr.abnorm <- c("108", "545", "564", "715", "B36", "HV80", "N7", "S23", "SE15285", "BJ17183", "AL9890", "506", "769", "887", "841", "903", "961", "1060", "1066", "TL14516", "14197")

#----
# ADD ATTRIBUTES TO CLINICAL TABLE
#----

c$jak2 <- ifelse(c$patient_id %in% patients.nonrelapsing.nomutdata, NA, c$patient_id %in% m[m$gene=="JAK2", "patient"])
c$jak2.dia <- ifelse(c$patient_id %in% patients.nonrelapsing.nomutdata, NA, c$patient_id %in% m[m$gene=="JAK2" & m$sample=="rem_dia", "patient"])
c$jak2.relapsing.dia <- ifelse(c$patient_id %in% patients.relapsing, c$jak2.dia, NA)
c$jak2.relapsing.rel <- ifelse(c$patient_id %in% patients.relapsing.matched, (c$patient_id %in% m[m$gene=="JAK2" & (m$sample=="rem_rel" | m$sample=="rem_rel2"), "patient"]), NA)
c$jak2.relapsing <- c$jak2.relapsing.dia | c$jak2.relapsing.rel 

c$ikzf1 <- c$patient_id %in% c(ikzf1.del.dia, ikzf1.del.rel)
c$ikzf1.dia <- c$patient_id %in% ikzf1.del.dia
c$ikzf1.relapsing.dia <- ifelse(c$patient_id %in% patients.relapsing, c$ikzf1.dia, NA)
c$ikzf1.relapsing.rel <- ifelse(c$patient_id %in% patients.relapsing.matched, c$patient_id %in% ikzf1.del.rel, NA)
c$ikzf1.relapsing <- c$ikzf1.relapsing.dia | c$ikzf1.relapsing.rel

c$cdkn2a <- c$patient_id %in% c(cdkn2a.del.dia, cdkn2a.del.rel)
c$cdkn2a.dia <- c$patient_id %in% cdkn2a.del.dia
c$cdkn2a.relapsing.dia <- ifelse(c$patient_id %in% patients.relapsing, c$cdkn2a.dia, NA)
c$cdkn2a.relapsing.rel <- ifelse(c$patient_id %in% patients.relapsing.matched, c$patient_id %in% cdkn2a.del.rel, NA)
c$cdkn2a.relapsing <- c$cdkn2a.relapsing.dia | c$cdkn2a.relapsing.rel

c$kras <- ifelse(c$patient_id %in% patients.nonrelapsing.nomutdata, NA, c$patient_id %in% m[m$gene=="KRAS", "patient"])
c$kras.dia <- ifelse(c$patient_id %in% patients.nonrelapsing.nomutdata, NA, c$patient_id %in% m[m$gene=="KRAS" & m$sample=="rem_dia", "patient"])
c$kras.relapsing.dia <- ifelse(c$patient_id %in% patients.relapsing, c$kras.dia, NA)
c$kras.relapsing.rel <- ifelse(c$patient_id %in% patients.relapsing.matched, (c$patient_id %in% m[m$gene=="KRAS" & (m$sample=="rem_rel" | m$sample=="rem_rel2"), "patient"]), NA)
c$kras.relapsing <- c$kras.relapsing.dia | c$kras.relapsing.rel 

lymphoid_devel.genes <- c("IKZF1", "IKZF2", "IKZF3", "PAX5", "EBF1", "ETV6")
c$lymphoid_devel <- c$patient_id %in% c(m[m$gene %in% lymphoid_devel.genes, "patient"], ikzf1.del.dia, ikzf1.del.rel, ikzf2.del.dia, ikzf2.del.rel, ikzf3.del.dia, ikzf3.del.rel, pax5.del.dia, pax5.del.rel, ebf1.del.dia, ebf1.del.rel, etv6.del.dia, etv6.del.rel) 
c$lymphoid_devel.dia <- c$patient_id %in% c(m[m$gene %in% lymphoid_devel.genes & m$sample=="rem_dia", "patient"], ikzf1.del.dia, ikzf2.del.dia, ikzf3.del.dia, pax5.del.dia, ebf1.del.dia, etv6.del.dia) 
c$lymphoid_devel.relapsing.dia <- ifelse(c$patient_id %in% patients.relapsing, c$lymphoid_devel.dia, NA)
c$lymphoid_devel.relapsing.rel <- ifelse(c$patient_id %in% patients.relapsing.matched, c$patient_id %in% c(m[m$gene %in% lymphoid_devel.genes & m$sample=="rem_rel", "patient"], ikzf1.del.rel, ikzf2.del.rel, ikzf3.del.rel, pax5.del.rel, ebf1.del.rel, etv6.del.rel), NA)
c$lymphoid_devel.relapsing <- c$lymphoid_devel.relapsing.dia | c$lymphoid_devel.relapsing.rel

jak_stat.genes <- c("JAK2", "JAK3", "JAK1", "IL7R", "SYK", "CRLF2")
c$jak_stat <- ifelse(c$patient_id %in% patients.nonrelapsing.nomutdata, NA, c$patient_id %in% m[m$gene %in% jak_stat.genes, "patient"]) 
c$jak_stat.dia <- ifelse(c$patient_id %in% patients.nonrelapsing.nomutdata, NA, c$patient_id %in% m[m$gene %in% jak_stat.genes & m$sample=="rem_dia", "patient"]) 
c$jak_stat.relapsing.dia <- ifelse(c$patient_id %in% patients.relapsing, c$jak_stat.dia, NA)
c$jak_stat.relapsing.rel <- ifelse(c$patient_id %in% patients.relapsing.matched, c$patient_id %in% m[m$gene %in% jak_stat.genes & m$sample=="rem_rel", "patient"], NA)
c$jak_stat.relapsing <- c$jak_stat.relapsing.dia | c$jak_stat.relapsing.rel

ras_pw.genes <- c("KRAS", "NRAS", "PTPN11", "FLT3")
c$ras_pw <- ifelse(c$patient_id %in% patients.nonrelapsing.nomutdata, NA, c$patient_id %in% m[m$gene %in% ras_pw.genes, "patient"])
c$ras_pw.dia <- ifelse(c$patient_id %in% patients.nonrelapsing.nomutdata, NA, c$patient_id %in% m[m$gene %in% ras_pw.genes & m$sample=="rem_dia", "patient"])
c$ras_pw.relapsing.dia <- ifelse(c$patient_id %in% patients.relapsing, c$ras_pw.dia, NA)
c$ras_pw.relapsing.rel <- ifelse(c$patient_id %in% patients.relapsing.matched, (c$patient_id %in% m[m$gene %in% ras_pw.genes & (m$sample=="rem_rel" | m$sample=="rem_rel2"), "patient"]), NA)
c$ras_pw.relapsing <- c$ras_pw.relapsing.dia | c$ras_pw.relapsing.rel 

c$signaling <- c$jak_stat | c$ras_pw
c$signaling.dia <- c$jak_stat.dia | c$ras_pw.dia
c$signaling.relapsing.dia <- c$jak_stat.relapsing.dia | c$ras_pw.relapsing.dia
c$signaling.relapsing.rel <- c$jak_stat.relapsing.rel | c$ras_pw.relapsing.rel
c$signaling.relapsing <- c$jak_stat.relapsing | c$ras_pw.relapsing

epigenetic.genes <- c("TRRAP", "SETD2", "CREBBP", "EP300", "MSH6")
c$epigenetic <- ifelse(c$patient_id %in% patients.nonrelapsing.nomutdata, NA, c$patient_id %in% c(m[m$gene %in% epigenetic.genes, "patient"], setd2.del.dia, setd2.del.rel)) 
c$epigenetic.dia <- ifelse(c$patient_id %in% patients.nonrelapsing.nomutdata, NA, c$patient_id %in% c(m[m$gene %in% epigenetic.genes & m$sample=="rem_dia", "patient"], setd2.del.dia, setd2.del.rel)) 
c$epigenetic.relapsing.dia <- ifelse(c$patient_id %in% patients.relapsing, c$epigenetic.dia, NA)
c$epigenetic.relapsing.rel <- ifelse(c$patient_id %in% patients.relapsing.matched, c$patient_id %in% c(m[m$gene %in% epigenetic.genes & m$sample=="rem_rel", "patient"], setd2.del.dia, setd2.del.rel), NA)
c$epigenetic.relapsing <- c$epigenetic.relapsing.dia | c$epigenetic.relapsing.rel

c$pax5 <- c$patient_id %in% c(pax5.del.dia, pax5.del.rel)
c$pax5[c$patient %in% c("GI8")] <- NA
c$pax5.dia <- c$patient_id %in% pax5.del.dia
c$pax5.dia[c$patient %in% c("GI8")] <- NA
c$pax5.relapsing.dia <- ifelse(c$patient_id %in% patients.relapsing, c$pax5.dia, NA)
c$pax5.relapsing.rel <- ifelse(c$patient_id %in% patients.relapsing.matched, c$patient_id %in% pax5.del.rel, NA)
c$pax5.relapsing.rel[c$patient %in% c("GI8")] <- NA
c$pax5.relapsing <- c$pax5.relapsing.dia | c$pax5.relapsing.rel

c$ikzf1.plus <- c$ikzf1 & (c$cdkn2a | c$pax5)
c$ikzf1.plus.dia <- c$ikzf1.dia & (c$cdkn2a.dia | c$pax5.dia)
c$ikzf1.plus.relapsing.dia <- c$ikzf1.relapsing.dia & (c$cdkn2a.relapsing.dia | c$pax5.relapsing.dia)
c$ikzf1.plus.relapsing.rel <- c$ikzf1.relapsing.rel & (c$cdkn2a.relapsing.rel | c$pax5.relapsing.rel)
c$ikzf1.plus.relapsing <- c$ikzf1.relapsing & (c$cdkn2a.relapsing | c$pax5.relapsing)

		
c$p2ry8.conserved <- ifelse(c$patient_id %in% patients.relapsing.matched, c$patient_id %in% p2ry8.cons, NA)

c$chr21.gain.not.iamp  <- c$patient_id %in% chr21.gain
c$chr21.gain.or.iamp <- c$chr21.gain.not.iamp | c$iAMP21
c$chr21.gain.or.iamp.somatic <- c$chr21.gain.or.iamp & !c$DS
#c$chr21.gain.or.iamp.somatic <- ifelse(c$DS, NA, c$chr21.gain.or.iamp)

c$sex.chr.abnorm <- c$patient_id %in% sex.chr.abnorm

write.table(c, file="~/p2ry8-crlf2/results/clinical/clinical_data.processed.tsv", col.names=T, row.names=F, sep="\t", quote=F)

#----
# ASSOCIATION TESTING
#----

pdf("~/p2ry8-crlf2/results/clinical/clinical-associations.significant.pdf")
tests <- test_pairwise_assoc(c, 
		sig.level=0.1, 
		exclude=c("patient_id", "source", "exome", "rel_protocol", "outcome_after.1st.rel", "BM.transplantation.date", "comment"),
		exclude.group=list(c("jak2", "jak2.dia", "jak2.relapsing.dia", "jak2.relapsing.rel", "jak2.relapsing", "jak_stat", "jak_stat.dia", "jak_stat.relapsing.dia", "jak_stat.relapsing.rel", "jak_stat.relapsing", "signaling", "signaling.dia", "signaling.relapsing.dia", "signaling.relapsing.rel", "signaling.relapsing"),
						   c("kras", "kras.dia", "kras.relapsing.dia", "kras.relapsing.rel", "kras.relapsing", "ras_pw", "ras_pw.dia", "ras_pw.relapsing.dia", "ras_pw.relapsing.rel", "ras_pw.relapsing", "signaling", "signaling.dia", "signaling.relapsing.dia", "signaling.relapsing.rel", "signaling.relapsing"),
						   c("ikzf1", "ikzf1.dia", "ikzf1.relapsing.dia", "ikzf1.relapsing.rel", "ikzf1.relapsing", "lymphoid_devel", "lymphoid_devel.dia", "lymphoid_devel.relapsing.dia", "lymphoid_devel.relapsing.rel", "lymphoid_devel.relapsing"),
						   c("cdkn2a", "cdkn2a.dia", "cdkn2a.relapsing.dia", "cdkn2a.relapsing.rel", "cdkn2a.relapsing"),
						   c("epigenetic", "epigenetic.dia", "epigenetic.relapsing.dia", "epigenetic.relapsing.rel", "epigenetic.relapsing"),
						   c("pax5", "pax5.dia", "pax5.relapsing.dia", "pax5.relapsing.rel", "pax5.relapsing", "lymphoid_devel", "lymphoid_devel.dia", "lymphoid_devel.relapsing.dia", "lymphoid_devel.relapsing.rel", "lymphoid_devel.relapsing"),
						   c("ikzf1.plus", "ikzf1.plus.dia", "ikzf1.plus.relapsing.dia", "ikzf1.plus.relapsing.rel", "ikzf1.plus.relapsing", "ikzf1", "ikzf1.dia", "ikzf1.relapsing.dia", "ikzf1.relapsing.rel", "ikzf1.relapsing", "cdkn2a", "cdkn2a.dia", "cdkn2a.relapsing.dia", "cdkn2a.relapsing.rel", "cdkn2a.relapsing", "pax5", "pax5.dia", "pax5.relapsing.dia", "pax5.relapsing.rel", "pax5.relapsing", "lymphoid_devel", "lymphoid_devel.dia", "lymphoid_devel.relapsing.dia", "lymphoid_devel.relapsing.rel", "lymphoid_devel.relapsing"),
						   c("chr21.gain.not.iamp", "chr21.gain.or.iamp", "iAMP21", "chr21.gain.or.iamp.somatic"),
						   c("chr21.gain.not.iamp", "chr21.gain.or.iamp", "chr21.gain.or.iamp.somatic", "DS"),
						   c("kras.relapsing.rel", "epigenetic"))
)
dev.off()


