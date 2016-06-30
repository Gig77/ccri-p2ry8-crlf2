t <- as.table(matrix(c(15, 8, 17, 1), nrow=2, dimnames = list('PrimaryLesionDia' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(13, 10, 13, 5), nrow=2, dimnames = list('SexChrAbnDia' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(16, 7, 13, 5), nrow=2, dimnames = list('IKZF1Dia' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(14, 9, 15, 3), nrow=2, dimnames = list('PAX5Dia' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(8, 15, 9, 9), nrow=2, dimnames = list('LymphoidDevelDia' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(13, 10, 13, 5), nrow=2, dimnames = list('JAK2Dia' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(9, 14, 11, 7), nrow=2, dimnames = list('JAK/STATDia' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(15, 8, 17, 1), nrow=2, dimnames = list('KRASDia' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(12, 11, 17, 1), nrow=2, dimnames = list('RTK/RasDia' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(14, 9, 14, 4), nrow=2, dimnames = list('EpigeneticRegDia' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(14, 9, 13, 5), nrow=2, dimnames = list('CDKN2ADia' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(12, 11, 13, 5), nrow=2, dimnames = list('OG/TSDia' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(8, 2, 5, 4), nrow=2, dimnames = list('JAK2Rel' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(7, 3, 4, 5), nrow=2, dimnames = list('JAK/STATRel' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(4, 6, 7, 2), nrow=2, dimnames = list('KRASRel' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(1, 9, 5, 4), nrow=2, dimnames = list('RTK/RASRel' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(7, 3, 9, 0), nrow=2, dimnames = list('CREBBPRel' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(7, 3, 9, 0), nrow=2, dimnames = list('KMT2DRel' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(5, 5, 8, 1), nrow=2, dimnames = list('EpigeneticRegRel' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(3, 7, 9, 1), nrow=2, dimnames = list('CDKN2ARel' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)

t <- as.table(matrix(c(1, 9, 7, 2), nrow=2, dimnames = list('OG/TSRel' = c("-", "+"), 'DS' = c('-', '+')))) ; t
sprintf("p=%.2g", fisher.test(t)$p.value)
