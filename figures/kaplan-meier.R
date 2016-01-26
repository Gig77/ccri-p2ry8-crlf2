options(warn=1)
library(survival)
library(cmprsk)

rm(list=ls())
c <- read.delim("/mnt/projects/p2ry8-crlf2/data/clinical_data_for_kaplan_meier.txt")
c$had_event <- c$event != "CCR"
c$is_dead <- !is.na(c$date.of.death) & c$date.of.death != ""
c$IKZF1_status_dia <- as.factor(c$IKZF1_status_dia)

# ==========================================================================================================

pdf("/mnt/projects/p2ry8-crlf2/results/figures/kaplan-maier.pdf", width=10, height=13.3)

par(mfrow=c(4,3),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))

# --- OS ---

form <- Surv(time=on_study_time_OS, is_dead)~IKZF1_status_dia
fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue"), lty=c(1, 1), xlab="months", ylab="pOS", conf.int=F, xlim=c(0,125), yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]-0.08, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0, cex=1.3)
legend("bottomright", c(sprintf("IKZF1 wt (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("IKZF1 mut (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("blue", "red"), lty=c(1, 1, 1), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

# --- EFS ---

rm(form) ; form <- Surv(time=on_study_time_EFS, had_event)~IKZF1_status_dia
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue"), lty=c(1, 1), xlab="months", ylab="pEFS", conf.int=F, xlim=c(0,125), yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0, cex=1.3)
legend("bottomright", c(sprintf("IKZF1 wt (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("IKZF1 mut (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("blue", "red"), lty=c(1, 1, 1), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

# cumulative incidence
# not done b/c there was no competing risk (all first events after diagnosis are relapses, no deaths)

#fit <- cuminc(c$second_rem_months, c$second_event_after_first_relapse, c$mrd_risk_rel, cencode="none")
#plot(fit[1:2], xlab="months", ylab="CIR", curvlab=c("MRD HR", "MRD SR"), col=c("red", "blue"), wh=c(-100,-100), yaxt='n', lty=c(1, 1), xlim=c(0,125), cex.axis=1.3, cex.lab=1.5)
#axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
#text(60, timepoints(fit, 60)$est["HR 2nd rel",]+0.03, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["HR 2nd rel",], timepoints(fit, 60)$var["HR 2nd rel",]), col="red", adj=0, cex=1.3)
#text(60, timepoints(fit, 60)$est["SR 2nd rel",]+0.03, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["SR 2nd rel",], timepoints(fit, 60)$var["SR 2nd rel",]), col="blue", adj=0, cex=1.3)
#box()
#legend("topright", c(sprintf("MRD HR (%d/%d)", sum(complete.cases(c[c$mrd_risk_rel=="HR" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$mrd_risk_rel=="HR",c("second_rem_months", "second_event_after_first_relapse")]))),
#                     sprintf("MRD SR (%d/%d)", sum(complete.cases(c[c$mrd_risk_rel=="SR" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$mrd_risk_rel=="SR",c("second_rem_months", "second_event_after_first_relapse")])))), 
#       lwd=c(1,1), col=c("red", "blue"), lty=c(1, 1, 1), box.lwd=0.5)
#text(1, 1, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", fit$Tests["2nd rel", "pv"])))), adj=0, cex=1.3)

dev.off()
