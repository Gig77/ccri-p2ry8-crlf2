# read data
d <- read.delim("/mnt/projects/p2ry8-crlf2/data/WES_CNADetection_Performance.txt", skip = 1, stringsAsFactors = F, na.strings = "")
names(d) <- c("Sample", "Cohort", "PAR1.SNP", "PAR1.PCR", "PAR1.WES", "IKZF1.SNP", "IKZF1.MLPA", "IKZF1.WES", "CDKN2A.SNP", "CDKN2A.MLPA", "CDKN2A.WES", "PAX5.SNP", "PAX5.FISH", "PAX5.WES", "Comment")

# setup result dataframe
perf <- data.frame("Gene" = character(0), "Metric" = character(0), "Percent" = numeric(0))

# PAR1 performance
par1 <- d[,c("Sample", "PAR1.SNP", "PAR1.PCR", "PAR1.WES")]
par1 <- par1[!is.na(par1$PAR1.WES) & (!is.na(par1$PAR1.SNP) | !is.na(par1$PAR1.PCR)),]
par1 <- par1[par1$PAR1.PCR != "SC",]
par1.tp <- sum(par1$PAR1.WES %in% c("del") & (par1$PAR1.SNP %in% c("del") | par1$PAR1.PCR %in% "del"))
par1.fp <- sum(par1$PAR1.WES %in% c("del") & ((is.na(par1$PAR1.SNP) | par1$PAR1.SNP %in% c("wt")) & (is.na(par1$PAR1.PCR) | par1$PAR1.PCR %in% "wt")))
par1.tn <- sum(par1$PAR1.WES %in% c("wt") &  ((is.na(par1$PAR1.SNP) | par1$PAR1.SNP %in% c("wt")) & (is.na(par1$PAR1.PCR) | par1$PAR1.PCR %in% "wt")))
par1.fn <- sum(par1$PAR1.WES %in% c("wt") &  (par1$PAR1.SNP %in% c("del") | par1$PAR1.PCR %in% "del"))
perf <- rbind(perf, setNames(data.frame("PAR1", "Sensitivity", par1.tp/(par1.tp+par1.fn), stringsAsFactors = F), names(perf)))
perf <- rbind(perf, setNames(data.frame("PAR1", "Specificity", par1.tn/(par1.tn+par1.fp), stringsAsFactors = F), names(perf)))
perf <- rbind(perf, setNames(data.frame("PAR1", "Accuracy", (par1.tp+par1.tn)/(par1.tp+par1.tn+par1.fp+par1.fn), stringsAsFactors = F), names(perf)))

# IKZF1 performance
ikzf1 <- d[,c("Sample", "IKZF1.SNP", "IKZF1.MLPA", "IKZF1.WES")]
ikzf1 <- ikzf1[!is.na(ikzf1$IKZF1.WES) & (!is.na(ikzf1$IKZF1.SNP) | !is.na(ikzf1$IKZF1.MLPA)),]
ikzf1.tp <- sum(ikzf1$IKZF1.WES %in% c("del") & (ikzf1$IKZF1.SNP %in% c("del") | ikzf1$IKZF1.MLPA %in% "del"))
ikzf1.fp <- sum(ikzf1$IKZF1.WES %in% c("del") & ((is.na(ikzf1$IKZF1.SNP) | ikzf1$IKZF1.SNP %in% c("wt")) & (is.na(ikzf1$IKZF1.MLPA) | ikzf1$IKZF1.MLPA %in% "wt")))
ikzf1.tn <- sum(ikzf1$IKZF1.WES %in% c("wt") &  ((is.na(ikzf1$IKZF1.SNP) | ikzf1$IKZF1.SNP %in% c("wt")) & (is.na(ikzf1$IKZF1.MLPA) | ikzf1$IKZF1.MLPA %in% "wt")))
ikzf1.fn <- sum(ikzf1$IKZF1.WES %in% c("wt") &  (ikzf1$IKZF1.SNP %in% c("del") | ikzf1$IKZF1.MLPA %in% "del"))
perf <- rbind(perf, setNames(data.frame("IKZF1", "Sensitivity", ikzf1.tp/(ikzf1.tp+ikzf1.fn), stringsAsFactors = F), names(perf)))
perf <- rbind(perf, setNames(data.frame("IKZF1", "Specificity", ikzf1.tn/(ikzf1.tn+ikzf1.fp), stringsAsFactors = F), names(perf)))
perf <- rbind(perf, setNames(data.frame("IKZF1", "Accuracy", (ikzf1.tp+ikzf1.tn)/(ikzf1.tp+ikzf1.tn+ikzf1.fp+ikzf1.fn), stringsAsFactors = F), names(perf)))

# CDKN2A performance
cdkn2a <- d[,c("Sample", "CDKN2A.SNP", "CDKN2A.MLPA", "CDKN2A.WES")]
cdkn2a <- cdkn2a[!is.na(cdkn2a$CDKN2A.WES) & (!is.na(cdkn2a$CDKN2A.SNP) | !is.na(cdkn2a$CDKN2A.MLPA)),]
cdkn2a.tp <- sum(cdkn2a$CDKN2A.WES %in% c("del") & (cdkn2a$CDKN2A.SNP %in% c("del") | cdkn2a$CDKN2A.MLPA %in% "del"))
cdkn2a.fp <- sum(cdkn2a$CDKN2A.WES %in% c("del") & ((is.na(cdkn2a$CDKN2A.SNP) | cdkn2a$CDKN2A.SNP %in% c("wt")) & (is.na(cdkn2a$CDKN2A.MLPA) | cdkn2a$CDKN2A.MLPA %in% "wt")))
cdkn2a.tn <- sum(cdkn2a$CDKN2A.WES %in% c("wt") &  ((is.na(cdkn2a$CDKN2A.SNP) | cdkn2a$CDKN2A.SNP %in% c("wt")) & (is.na(cdkn2a$CDKN2A.MLPA) | cdkn2a$CDKN2A.MLPA %in% "wt")))
cdkn2a.fn <- sum(cdkn2a$CDKN2A.WES %in% c("wt") &  (cdkn2a$CDKN2A.SNP %in% c("del") | cdkn2a$CDKN2A.MLPA %in% "del"))
perf <- rbind(perf, setNames(data.frame("CDKN2A", "Sensitivity", cdkn2a.tp/(cdkn2a.tp+cdkn2a.fn), stringsAsFactors = F), names(perf)))
perf <- rbind(perf, setNames(data.frame("CDKN2A", "Specificity", cdkn2a.tn/(cdkn2a.tn+cdkn2a.fp), stringsAsFactors = F), names(perf)))
perf <- rbind(perf, setNames(data.frame("CDKN2A", "Accuracy", (cdkn2a.tp+cdkn2a.tn)/(cdkn2a.tp+cdkn2a.tn+cdkn2a.fp+cdkn2a.fn), stringsAsFactors = F), names(perf)))

# PAX5 performance
pax5 <- d[,c("Sample", "PAX5.SNP", "PAX5.FISH", "PAX5.WES")]
pax5 <- pax5[!is.na(pax5$PAX5.WES) & (!is.na(pax5$PAX5.SNP) | !is.na(pax5$PAX5.FISH)),]
pax5.tp <- sum(pax5$PAX5.WES %in% c("del") & (pax5$PAX5.SNP %in% c("del") | pax5$PAX5.FISH %in% "del"))
pax5.fp <- sum(pax5$PAX5.WES %in% c("del") & ((is.na(pax5$PAX5.SNP) | pax5$PAX5.SNP %in% c("wt")) & (is.na(pax5$PAX5.FISH) | pax5$PAX5.FISH %in% "wt")))
pax5.tn <- sum(pax5$PAX5.WES %in% c("wt") &  ((is.na(pax5$PAX5.SNP) | pax5$PAX5.SNP %in% c("wt")) & (is.na(pax5$PAX5.FISH) | pax5$PAX5.FISH %in% "wt")))
pax5.fn <- sum(pax5$PAX5.WES %in% c("wt") &  (pax5$PAX5.SNP %in% c("del") | pax5$PAX5.FISH %in% "del"))
perf <- rbind(perf, setNames(data.frame("PAX5", "Sensitivity", pax5.tp/(pax5.tp+pax5.fn), stringsAsFactors = F), names(perf)))
perf <- rbind(perf, setNames(data.frame("PAX5", "Specificity", pax5.tn/(pax5.tn+pax5.fp), stringsAsFactors = F), names(perf)))
perf <- rbind(perf, setNames(data.frame("PAX5", "Accuracy", (pax5.tp+pax5.tn)/(pax5.tp+pax5.tn+pax5.fp+pax5.fn), stringsAsFactors = F), names(perf)))

# overall performance
overall.tp <- par1.tp + ikzf1.tp + cdkn2a.tp + pax5.tp
overall.fp <- par1.fp + ikzf1.fp + cdkn2a.fp + pax5.fp
overall.tn <- par1.tn + ikzf1.tn + cdkn2a.tn + pax5.tn
overall.fn <- par1.fn + ikzf1.fn + cdkn2a.fn + pax5.fn
perf <- rbind(perf, setNames(data.frame("Overall", "Sensitivity", overall.tp/(overall.tp+overall.fn), stringsAsFactors = F), names(perf)))
perf <- rbind(perf, setNames(data.frame("Overall", "Specificity", overall.tn/(overall.tn+overall.fp), stringsAsFactors = F), names(perf)))
perf <- rbind(perf, setNames(data.frame("Overall", "Accuracy", (overall.tp+overall.tn)/(overall.tp+overall.tn+overall.fp+overall.fn), stringsAsFactors = F), names(perf)))

# function for grouped bar plot using base graphics
# copied from https://github.com/mrxiaohe/R_Functions/blob/master/functions/bar
bar <- function(dv, factors, dataframe, percentage=FALSE, errbar=!percentage, half.errbar=TRUE, conf.level=.95, 
                xlab=NULL, ylab=NULL, main=NULL, names.arg=NULL, bar.col="black", whisker=.015,args.errbar=NULL,
                legend=TRUE, legend.text=NULL, args.legend=NULL,legend.border=FALSE, box=TRUE, args.yaxis=NULL, 
                mar=c(5,4,3,2),...){
  axes=!percentage
  dv.name<-substitute(dv)
  if(length(dv.name)>1) stop("'dv' only takes one variable")
  dv.name<-as.character(dv.name)
  dv<-dataframe[[dv.name]]
  fnames<-substitute(factors)
  if(length(fnames)==1){
    factors<-as.character(fnames)
    nf<-1
  }else{
    factors<-as.character(fnames[-1L])
    nf<-length(factors)
  }
  if(nf>2) stop("This function accepts no more than 2 factors \n",
                "\t-i.e., it only plots one-way or two-way designs.")
  if(percentage & errbar){
    warning("percentage=TRUE; error bars were not plotted")
    errbar<-FALSE
  }
  if(!percentage) xbars<-tapply(dv, dataframe[,factors], mean, na.rm=TRUE)
  else {
    xbars<-tapply(dv, list(interaction(dataframe[,factors], lex.order=TRUE)), mean, na.rm=TRUE)
    if(sum(na.omit(dv)!=0&na.omit(dv)!=1)>0) 
      stop("Data points in 'dv' need to be 0 or 1 in order to set 'percentage' to TRUE")
    xbars<-rbind(xbars, 1-xbars)*100
  }
  if(errbar){
    se<-tapply(dv, dataframe[,factors], sd, na.rm=TRUE)/sqrt(tapply(dv, dataframe[,factors], length))
    conf.level=1-(1-conf.level)/2
    lo.bar<-xbars-se*qnorm(conf.level)
    hi.bar<-xbars+se*qnorm(conf.level)	
  }
  extras<-list(...)
  if(legend & !percentage){
    if(is.null(legend.text))
      legend.text<-sort(unique(dataframe[[factors[1]]]))
    args.legend.temp<-list(x="topright", bty=if(!legend.border)"n" else "o",
                           inset=c(0,0))
    if(is.list(args.legend))
      args.legend<-modifyList(args.legend.temp, args.legend)
    else 
      args.legend<-args.legend.temp
  } else if(legend & percentage){
    if(is.null(legend.text)) 
      legend.text<-c("1", "0")
    args.legend.temp<-list(x="topright", bty=if(!legend.border)"n" else "o",
                           inset=c(0,0))
    if(is.list(args.legend))
      args.legend<-modifyList(args.legend.temp, args.legend)
    else 
      args.legend<-args.legend.temp
  } else if(!legend){
    args.legend<-NULL
    legend.text<-NULL
  }
  if(errbar && legend && !percentage) ymax<-max(hi.bar)+max(hi.bar)/20
  else if(errbar && legend && percentage) ymax<-115
  else if(errbar && !legend) ymax <- max(xbars)
  else if(!errbar && legend && percentage) ymax<-110	
  else if(!errbar) ymax<-max(xbars) + max(xbars)/20
  if(!percentage){
    args.barplot<-list(beside=TRUE, height=xbars, ylim=c(0, ymax), main=main, names.arg=names.arg,
                       col=hcl(h=seq(0,270, 270/(length(unique(dataframe[[factors[1]]]))))[-length(unique(dataframe[[factors[1]]]))]),
                       legend.text=legend.text, args.legend=args.legend, xpd=TRUE,
                       xlab=if(is.null(xlab)) factors[length(factors)] else xlab,
                       ylab=if(is.null(ylab)) dv.name else ylab, axes=axes)
  }else{
    args.barplot<-list(beside=TRUE, height=xbars, ylim=c(0, ymax),  main=main, names.arg=names.arg,
                       col=hcl(h=seq(0,270, 270/(length(unique(dataframe[[factors[1]]]))))[-length(unique(dataframe[[factors[1]]]))]),
                       legend.text=legend.text, args.legend=args.legend, xpd=TRUE,
                       xlab=if(is.null(xlab)) " "[length(factors)] else xlab,
                       ylab=if(is.null(ylab)) "percentage" else ylab, axes=axes)		
  }
  args.barplot<-modifyList(args.barplot, extras)
  errbars = function(xvals, cilo, cihi, whisker, nc, args.errbar = NULL, half.errbar=TRUE) {
    if(half.errbar){
      cilo<-(cihi+cilo)/2
    }
    fixedArgs.bar = list(matlines, x=list(xvals), 
                         y=lapply(split(as.data.frame(t(do.call("rbind", 
                                                                list(cihi, cilo)))),1:nc),matrix, 
                                  nrow=2, byrow=T))
    allArgs.bar = c(fixedArgs.bar, args.errbar)
    whisker.len = whisker*(par("usr")[2] - par("usr")[1])/2
    whiskers = rbind((xvals - whisker.len)[1,],
                     (xvals + whisker.len)[1,])
    fixedArgs.lo = list(matlines, x=list(whiskers), 	
                        y=lapply(split(as.data.frame(t(do.call("rbind", 
                                                               list(cilo, cilo)))), 1:nc), matrix, nrow=2, byrow=T))
    allArgs.bar.lo = c(fixedArgs.lo, args.errbar)
    fixedArgs.hi = list(matlines, x=list(whiskers), 
                        y=lapply(split(as.data.frame(t(do.call("rbind", 
                                                               list(cihi, cihi)))), 1:nc), matrix, nrow=2, byrow=T))
    allArgs.bar.hi = c(fixedArgs.hi, args.errbar)  
    invisible(do.call(mapply, allArgs.bar))
    if(!half.errbar) invisible(do.call(mapply, allArgs.bar.lo))
    invisible(do.call(mapply, allArgs.bar.hi))
  }
  par(mar=mar)
  errloc<-as.vector(do.call(barplot, args.barplot))
  if(errbar){
    errloc<-rbind(errloc, errloc)
    lo.bar<-matrix(as.vector(lo.bar))
    hi.bar<-matrix(as.vector(hi.bar))
    args.errbar.temp<-list(col=bar.col, lty=1)
    args.errbar<-if(is.null(args.errbar)|!is.list(args.errbar)) 
      args.errbar.temp
    else if(is.list(args.errbar)) 
      modifyList(args.errbar.temp, args.errbar)
    errbars(errloc, cilo=lo.bar, cihi=hi.bar, nc=1, whisker=whisker, 
            args.errbar=args.errbar, half.errbar=half.errbar)
  }
  if(box) box()
  if(percentage){
    args.yaxis.temp<-list(at=seq(0,100, 20), las=1)
    args.yaxis<-if(!is.list(args.yaxis)) args.yaxis.temp else modifyList(args.yaxis.temp, args.yaxis)
    do.call(axis, c(side=2, args.yaxis))
  }
}

# barplot
perf$Gene <- factor(perf$Gene, levels=c("PAR1", "IKZF1", "CDKN2A", "PAX5", "Overall"))
perf$Metric <- factor(perf$Metric, levels=c("Sensitivity", "Specificity", "Accuracy"))
cols <- c("red", "blue", "black")

pdf("/mnt/projects/p2ry8-crlf2/results/figures/wes-cnv-performance.pdf", height=6.6)
par(xpd=T)
bar(Percent, factors = c(Metric, Gene), perf, errbar = FALSE, ylim=c(0, 1.15), col=cols, legend=F, xlab="", mar=c(10,4.1,4.1,2.1), mgp=c(3, 0.5, 0), tck=-0.02)
legend("top", legend=levels(perf$Metric), fill=cols, ncol=3, bty="n")
rect(0.2, -0.45, 20.8, -0.15)
y <- -0.2 ; text(-1.2, y, "TP", adj=c(0, 0.5)) ; text(2.5, y, par1.tp) ; text(6.5, y, ikzf1.tp) ; text(10.5, y, cdkn2a.tp) ; text(14.5, y, pax5.tp) ; text(18.5, y, overall.tp)
y <- -0.27 ; text(-1.2, y, "FP", adj=c(0, 0.5)) ; text(2.5, y, par1.fp) ; text(6.5, y, ikzf1.fp) ; text(10.5, y, cdkn2a.fp) ; text(14.5, y, pax5.fp) ; text(18.5, y, overall.fp)
y <- -0.34 ; text(-1.2, y, "TN", adj=c(0, 0.5)) ; text(2.5, y, par1.tn) ; text(6.5, y, ikzf1.tn) ; text(10.5, y, cdkn2a.tn) ; text(14.5, y, pax5.tn) ; text(18.5, y, overall.tn)
y <- -0.41 ; text(-1.2, y, "FN", adj=c(0, 0.5)) ; text(2.5, y, par1.fn) ; text(6.5, y, ikzf1.fn) ; text(10.5, y, cdkn2a.fn) ; text(14.5, y, pax5.fn) ; text(18.5, y, overall.fn)
dev.off()

