
# DS

t <- as.table(matrix(c(9,10,9,13), nrow = 2, dimnames = list('DS'=c("yes", "no"), 'Cohort'=c("Non-relapsing", "Relapsing")))) ; t
fisher.test(t, workspace = 2e9)$p.value

t <- as.table(matrix(c(18,23,9,11), nrow = 2, dimnames = list('DS'=c("yes", "no"), 'Study'=c("This", "BFMAT")))) ; t
fisher.test(t, workspace = 2e9)$p.value

# age at diagnosis
age.nonrel <- c(3.8, 3.5, 3.0, 1.8, 2.4, 4.7, 4.6, 6.6, 5.7, 3.7, 3.9, 2.9, 3.2, 3.8, 5.6, 3.7, 5.8, 3.3, 4.4)
age.rel <- c(8.7, 4.1, 2.7, 4.6, 9.4, 9.6, 4.1, 2.0, 6.1, 1.5, 8.1, 4.6, 2.6, 3.0, 3.8, 8.9, 7.9, 2.6, 4.3, 1.8, 3.3, 7.6)
age.ctrl <- c(3.9, 8.5, 6.8, 13.6, 2.50, 3.50, 5.5, 2.5, 6.5, 2.2, 3.8, 18.2, 4.0, 4.0, 4.0, 6.7, 17, 14.1, 3.5, 3.4)

wilcox.test(age.nonrel, age.rel, alternative="two.sided", paired=FALSE)
wilcox.test(c(age.nonrel, age.rel), age.ctrl, alternative="two.sided", paired=FALSE)

# sex

t <- as.table(matrix(c(6,13,12,10), nrow = 2, dimnames = list('Sex'=c("Female", "Male"), 'Cohort'=c("Non-relapsing", "Relapsing")))) ; t
fisher.test(t, workspace = 2e9)$p.value

t <- as.table(matrix(c(18,23,7,13), nrow = 2, dimnames = list('Sex'=c("Female", "Male"), 'Study'=c("This", "BFMAT")))) ; t
fisher.test(t, workspace = 2e9)$p.value

# WCC

wcc.nonrel <- c(2, 86.3, 2.3, 30, 3.8, 7.8, 3.5, 5.5, 2.5, 114000, 3500, 6300, 1.7, 2.4, 5.2, 4.5, 19.4, 11.4, 2.3)
wcc.rel <- c(3.6, 69.2, 83.7, 5.6, 67.8, 55.2, 95, 11.8, 19.4, 98.5, 75.7, 21.5, 14.5, 2.7, 27.2, 6.2, 61.4, 9.3, 452, 300, 158, 10.9)
wcc.ctrl <- c(0.5, 0.8, 1.3, 1.5, 1.54, 3, 3.23, 3.7, 5.18, 7.8, 10.99, 14.25, 19.75, 27.09, 48.4, 189.8, 222)

wilcox.test(wcc.nonrel, wcc.rel, alternative="two.sided", paired=FALSE)
wilcox.test(c(wcc.nonrel, wcc.rel), wcc.ctrl, alternative="two.sided", paired=FALSE)

#t <- as.table(matrix(c(2,17,11,11), nrow = 2, dimnames = list('WCC'=c(">50", "<50"), 'Cohort'=c("Non-relapsing", "Relapsing")))) ; t
#fisher.test(t, workspace = 2e9)$p.value

#t <- as.table(matrix(c(13,28,2,18), nrow = 2, dimnames = list('WCC'=c(">50", "<50"), 'Study'=c("This", "BFMAT")))) ; t
#fisher.test(t, workspace = 2e9)$p.value

# MRD RG

t <- as.table(matrix(c(9,7,2,1,5,11,4,2), nrow = 4, dimnames = list('MRD RG'=c("SR", "IR", "HR", "na"), 'Cohort'=c("Non-relapsing", "Relapsing")))) ; t
fisher.test(t, workspace = 2e9)$p.value

t <- as.table(matrix(c(14,18,6,3,3,13,2,2), nrow = 4, dimnames = list('MRD RG'=c("SR", "IR", "HR", "na"), 'Study'=c("This", "BFMAT")))) ; t
fisher.test(t, workspace = 2e9)$p.value

# Treatment arm

t <- as.table(matrix(c(9,8,2,5,13,4), nrow = 3, dimnames = list('Treatment arm'=c("SR", "IR", "HR"), 'Cohort'=c("Non-relapsing", "Relapsing")))) ; t
fisher.test(t, workspace = 2e9)$p.value

t <- as.table(matrix(c(14,21,6,3,15,2), nrow = 3, dimnames = list('Treatment arm'=c("SR", "IR", "HR"), 'Study'=c("This", "BFMAT")))) ; t
fisher.test(t, workspace = 2e9)$p.value
