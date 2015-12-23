load("/mnt/projects/p2ry8-crlf2/results/exomeCopy/counts.bg.RData")
counts <- counts[seqnames(counts) %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"),]
d <- counts@values@unlistData@listData

# compute smoothed GC
d$GC.smoothed <- NA
for (c in unique(seqnames(counts))) {
  d$GC.smoothed[seqnames(counts)==c] <- caTools::runmean(d$GC[seqnames(counts)==c], 1000)
}
names(d) <- paste0("d", names(d))
fit <- lm(log(d887D+0.1, 2) ~ dGC.sq + d948C, data=d)
summary(fit)
counts$'887D.adj' <- fitted(fit)

# plot chr 3

chr <- "3" ; s.tum <- "GI8R" ; s.norm <- "GI8C" ; chr.excl <- c("11", "21", "X")

counts.chr <- counts[seqnames(counts)==chr,]
x <- 1:nrow(counts.chr)

par(mfrow=c(2,1))
t <- log((counts.chr[[s.tum]]+0.1), 2)
n <- log((counts.chr[[s.norm]]+0.1), 2)
ratio <- t-n
ratio.loess <- loess(ratio ~ x, span=0.1)
ratio.smooth <- predict(ratio.loess, x)
plot(1:length(t), ratio, ylim=c(-2, 2), main=paste(s.tum, "chr", chr), cex=0.2)
lines(x, ratio.smooth, col="orange", lwd=2)

#t <- counts.chr[["887D.adj"]]
#ratio <- t-n
ratio.loess <- loess(ratio.adjusted ~ x, span=0.1)
ratio.smooth <- predict(ratio.loess, x)
plot(1:length(t), ratio.adjusted, ylim=c(-2, 2), main=paste(s.tum, "chr", chr, "(GC adjusted)"))
lines(x, ratio.smooth, col="orange", lwd=2)
plot(x, d$dGC.smoothed[seqnames(counts)=="3"])

par(mfrow=c(2,1))
plot(1:length(t), ratio, ylim=c(-2, 2), main=paste(s.tum, "chr", chr), cex=0.2)
plot(1:length(t), ratio+4*d$dGC.sq[seqnames(counts)==chr], ylim=c(-2, 2), main=paste(s.tum, "chr", chr, "adjusted"), cex=0.2)

counts.nonzero <- counts[counts[[s.norm]]>=10,]
ratio.genome <- log((counts.nonzero[[s.tum]]+1) / (counts.nonzero[[s.norm]]+1), 2)
ratio.genome <- ratio.genome-mean(ratio.genome)
counts.pruned <- counts.nonzero
#counts.pruned <- counts.nonzero[ratio.genome >= -2 & ratio.genome <= 2,]
ratio.genome <- log((counts.pruned[[s.tum]]+1) / (counts.pruned[[s.norm]]+1), 2)
ratio.genome <- ratio.genome-mean(ratio.genome)

gc.centered <- counts.pruned$GC.sq-mean(counts.pruned$GC.sq)
gc.scaled <- scale(counts.nonzero$GC.smoothed^2)
fit <- lm(ratio.genome ~ gc.centered)
summary(fit)
coeffs <- coef(fit)
chr <- "1"
par(mfrow=c(2,1))
x <- (1:length(ratio.genome))[seqnames(counts.nonzero)==chr]
y <- ratio.genome[seqnames(counts.nonzero)==chr]
plot(x, y, ylim=c(-2, 2), main=paste(s.tum), cex=0.3, col=rgb(0, 0, 0, 1))
points(x, predict(loess(y ~ x, span=0.1)), type="l", col="orange", lwd=8)
points(x, predict(loess(y ~ x, span=0.01)), type="l", col="blue", lwd=3)
#plot(x, (ratio.genome-gc.scaled/1.2)[seqnames(counts.nonzero)==chr], ylim=c(-2, 2), main=paste(s.tum), cex=0.3, col=rgb(0, 0, 0, 0.3))
#plot(x, (ratio.genome-predict(fit.loess))[seqnames(counts.nonzero)==chr], ylim=c(-2, 2), main=paste(s.tum), cex=0.3, col=rgb(0, 0, 0, 0.3))

# plot GC along chromosomes
x <- (1:length(ratio.genome))[seqnames(counts.nonzero)==chr]
#y <- scale(counts.nonzero$GC.sq)[seqnames(counts.nonzero)==chr] 
y <- counts.nonzero$GC[seqnames(counts.nonzero)==chr]
plot(x, y, main="GC", cex=0.3, col=rgb(0, 0, 0, 1), ylim=c(0.2, 0.8))
points(x, predict(loess(y ~ x, span=0.1)), type="l", col="orange", lwd=8)
points(x, predict(loess(y ~ x, span=0.01)), type="l", col="blue", lwd=3)

# raw counts
x <- (1:length(ratio.genome))[seqnames(counts.nonzero)==chr]
y <- log(counts.nonzero[[s.tum]]+1, 2)[seqnames(counts.nonzero)==chr]
plot(x, y, cex=0.5, col=rgb(0, 0, 0, 1))
plot(counts.nonzero$GC[seqnames(counts.nonzero)==chr], counts.nonzero[[s.tum]][seqnames(counts.nonzero)==chr], cex=0.5, col=rgb(0, 0, 0, 1), ylim=c(0, 500))

# scatter
s.tum <- "m248-841-dia"
s.norm <- "841C"
counts.nonzero <- counts[counts[[s.norm]]>=10 & counts$GC > 0.1 & counts$GC < 0.9,]
gc.smoothed <- rep(NA, nrow(counts.nonzero))
for (c in unique(seqnames(counts.nonzero))) {
  gc.smoothed[seqnames(counts.nonzero)==c] <- caTools::runmean(counts.nonzero$GC[seqnames(counts.nonzero)==c], 50)
}
counts.nonzero$GC.smoothed <- gc.smoothed

ratio.genome <- log((counts.nonzero[[s.tum]]+1) / (counts.nonzero[[s.norm]]+1), 2)
ratio.genome <- ratio.genome-mean(ratio.genome)
#gc.centered <- counts.nonzero$GC-mean(counts.nonzero$GC)
gc.centered <- counts.nonzero$GC
fit1 <- lm(ratio.genome ~ gc.centered)
summary(fit1)
fit2 <- lm(ratio.genome ~ gc.centered + I(gc.centered^2))
summary(fit2)
anova(fit1, fit2)
fit3 <- lm(ratio.genome ~ gc.centered + I(gc.centered^2) + I(gc.centered^3))
summary(fit3)
anova(fit2, fit3)
fit.loess <- loess(ratio.genome ~ counts.nonzero$GC.smoothed)
fitfun.pol1 <- function(x) fit1$coefficient[1] + fit1$coefficient[2]*x 
fitfun.pol2 <- function(x) fit2$coefficient[1] + fit2$coefficient[2]*x + fit2$coefficient[3]*x^2
fitfun.pol3 <- function(x) fit3$coefficient[1] + fit3$coefficient[2]*x + fit3$coefficient[3]*x^2  + fit3$coefficient[4]*x^3
plot(gc.centered, ratio.genome, cex=0.1, col=rgb(0,0,0,0.3), xlim=c(0, 1))
curve(fitfun.pol1, col="orange", lwd=3, add=T)
curve(fitfun.pol2, col="blue", lwd=3, add=T)
curve(fitfun.pol3, col="brown", lwd=3, add=T)
sampl <- sample(1:length(gc.centered), 1000) ; sampl <- sampl[order(gc.centered[sampl])]
points(gc.centered[sampl], predict(fit.loess)[sampl], type="l", col="green", lwd=5)

# fit on subsample
library(MASS)
sampl <- sample(1:length(gc.centered), 20000) ; sampl <- sampl[order(gc.centered[sampl])]
x <- gc.centered[sampl]
y <- ratio.genome[sampl]
fit.loess <- loess(y ~ x, family="symmetric")
fit2 <- rlm(y ~ x + I(x^2))
plot(x, y, cex=0.3, col=rgb(0,0,0,0.3))
points(gc.centered[sampl], predict(fit2), type="l", col="orange", lwd=5)
points(gc.centered[sampl], predict(fit.loess), type="l", col="orange", lwd=5)
summary(fit2)
