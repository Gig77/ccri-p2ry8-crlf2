library(Gviz)
library(GenomicRanges)

library(exomeCopy)

cychp = list(
  "737D" = "/mnt/projects/p2ry8-crlf2/data/snparray/A1311-15.CEL._2007-2323.na33.cyhd.cychp.txt",
  "737R" = "/mnt/projects/p2ry8-crlf2/data/snparray/2009-2490_A1311-16.CEL._na33.cyhd.cychp.txt",
  "948D" = "/mnt/projects/p2ry8-crlf2/data/snparray/2010-3831._na33.cyhd.cychp.txt",
  "HV80R" = "/mnt/projects/p2ry8-crlf2/data/snparray/St.Anna_HV80-R.NA33.cyhd.cychp.txt"
)

regions = data.frame(patient=character(0), name=character(0), chr=character(0), start=numeric(0), end=numeric(0), genes=character(0), stringsAsFactors = F)
regions = rbind(regions, setNames(data.frame("737D", "PAR1", "X", 539250, 2484490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("737R", "IKZF1", "7", 49925794, 50867611, "IKZF1", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("948D", "PAX5", "9", 36430065, 37406492, "PAX5", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("HV80R", "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("737D", "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("737R", "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))

pdf("/mnt/projects/p2ry8-crlf2/results/figures/exomeCopy.pdf", width=5.5/2.54, height=4/2.54)

cols <- c("darkblue", "blue", "darkgray", "red", "darkred", "darkred", "darkred")
par(pch=19)
for (i in 1:nrow(regions)) {
  region <- regions[i,]

  # get array data
  cmd <- paste0("cat ", cychp[[region$patient]], " | perl /mnt/projects/p2ry8-crlf2/scripts/cychp2csv.pl ", region$chr, ":", region$start, "-", region$end)
  print(cmd)
  array <- read.table(pipe(cmd), header=T, stringsAsFactors = F)
  array$Chromosome <- paste0("chr", array$Chromosome)
  lrr.gr <- makeGRangesFromDataFrame(array[array$Type=="LRR",c("Chromosome", "Start", "Value")], start.field = "Start", end.field = "Start", keep.extra.columns = T)
  cn.gr <- makeGRangesFromDataFrame(array[array$Type=="CN",c("Chromosome", "Start", "End", "Value")], start.field = "Start", end.field = "End", keep.extra.columns = T)
  
  # encode CN state into LRR
  o <- findOverlaps(lrr.gr, cn.gr)
  states <- cn.gr[o@subjectHits]$Value
  lrr.gr <- GRanges(seqnames=seqnames(lrr.gr), ranges=ranges(lrr.gr), 
              mcols=data.frame(
                Gain   = as.numeric(ifelse(states >  2, lrr.gr$Value, NA)),
                Loss   = as.numeric(ifelse(states <  2, lrr.gr$Value, NA)),
                Normal = as.numeric(ifelse(states == 2, lrr.gr$Value, NA))
              ))

  # get WES data
  load(paste0("/mnt/projects/p2ry8-crlf2/results/exomeCopy/", region$patient, ".exomeCopy.fit.RData"))
  states <- fit[[1]][[region$chr]]@fx.par$S
  path <- as.vector(fit[[1]][[region$chr]]@path)
  ratio <- fit[[1]][[region$chr]]@O.norm
  wes.gr <- GRanges(seqnames=paste0("chr", region$chr), 
                    ranges=fit[[1]][[region$chr]]@ranges[[1]], 
                    mcols=data.frame(
                      Normal=ifelse(states[path] == 2, ratio, NA),
                      Loss=ifelse(states[path] < 2, ratio, NA),
                      Gain=ifelse(states[path] > 2, ratio, NA)
                    ))

  # configure tracks
  ideogram.track <- IdeogramTrack(genome = "hg19", chromosome = paste0("chr", region$chr), showId=F, showTitle=T, 
                                  name=paste0("chr", region$chr), fontcolor.title="black", background.title="white", rotation.title=0, fontsize=40)
  axis.track <- GenomeAxisTrack(labelPos="above", name=paste0("chr", region$chr), lwd=0.5, cex=0.25, col="black", fontcolor="black", fontface=2)
  biomTrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = paste0("chr", region$chr), 
                                      start = region$start, end = region$end, name = "", showId=TRUE, 
                                      stacking="squish", collapseTranscripts = "longest",
                                      fontsize=7, fontcolor.group="black", background.title="white")
  biomTrack@range <- biomTrack@range[biomTrack@range$symbol %in% unlist(strsplit(region$genes, ","))]
  lrr.track <- DataTrack(lrr.gr, name = paste(region$patient, "Array"), type = c("p", "g"), groups=c("Gain", "Loss", "Normal"), col=c("red", "blue", "darkgray"), 
                         cex=0.2, lty.grid=2, v=0, ylim=c(-1.8, 1.8), alpha.title = 1, alpha=0.8, fontsize=40, col.title="black", col.axis="black")
  sep.track <- DataTrack(ylim=c(1,1), showTitle=F, showAxis=F, background.panel = "lightgray")
  wes.track <- DataTrack(wes.gr, name = paste(region$patient, "WES"), type = c("p", "g"), groups=c("Gain", "Loss", "Normal"), col=c("darkgray", "blue", "red"), 
                         cex=0.4, lty.grid=2, v=0, ylim=c(0.1,1.9), alpha.title = 1, alpha=0.8, fontsize=40, col.title="black", col.axis="black")
  
  #pdf("/mnt/projects/p2ry8-crlf2/results/figures/exomeCopy.pdf", width=5.5/2.54, height=4/2.54)
  plotTracks(
    trackList = c(ideogram.track, axis.track, biomTrack, lrr.track, sep.track, wes.track), 
    sizes     = c(          0.09,       0.195,      0.09,      0.31,     0.005,       0.31),
    from      = region$start, 
    to        = region$end,
    cex.title = 0.1, cex.axis=0.1
  )
  #dev.off()
  
  next
  
  plot(
    x=fit[[region$patient]][[region$chr]], 
    main=paste(region$patient, region$name), 
    cex=0.7, 
    cex.main=1, 
    cex.axis=1, 
    ylim=c(0, 4), 
    xlim=c(region$start, region$end), 
    xlab=paste("Chromosome", region$chr), 
    ylab="Copy number", 
    col=cols,
    show.legend=F,
    panel.first=abline(h=1:4, lty=2, col="#dddddd", lwd=1)
  )
}

dev.off()
