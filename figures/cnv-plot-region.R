library(foreach)
library(doMC)
library(Gviz)
library(GenomicRanges)
library(exomeCopy)

registerDoMC(40) # number of cores for parallel execution

cychp = list(
  "108D"  = "NA",
  "108R"  = "NA",
  "108R2" = "NA",
  "737D"  = "/mnt/projects/p2ry8-crlf2/data/snparray/A1311-15.CEL._2007-2323.na33.cyhd.cychp.txt",
  "737R"  = "/mnt/projects/p2ry8-crlf2/data/snparray/2009-2490_A1311-16.CEL._na33.cyhd.cychp.txt",
  "737R2" = "NA",
  "737R3" = "/mnt/projects/p2ry8-crlf2/data/snparray/2012-2407._na33.cyhd.cychp.txt",
  "715D"  = "/mnt/projects/p2ry8-crlf2/data/snparray/2007-0966_A1311-12.CEL.NA33.cyhd.cychp.txt",
  "715R"  = "/mnt/projects/p2ry8-crlf2/data/snparray/2008-2775_A1311-13.CEL..NA33.cyhd.cychp.txt",
  "715R3" = "/mnt/projects/p2ry8-crlf2/data/snparray/715_RR._na33.cyhd.ND.cychp.txt",
  "92D"   = "/mnt/projects/p2ry8-crlf2/data/snparray/94-1115_92_D._na33.cyhd.cychp.txt",
  "92R"   = "/mnt/projects/p2ry8-crlf2/data/snparray/98-0366_92_R._na33.cyhd.cychp.txt",
  "839D"  = "/mnt/projects/p2ry8-crlf2/data/snparray/2009-1595_839_D._na33.cyhd.cychp.txt",
  "839R"  = "/mnt/projects/p2ry8-crlf2/data/snparray/2013-0667_839_R._na33.cyhd.cychp.txt",
  "948D"  = "/mnt/projects/p2ry8-crlf2/data/snparray/2010-3831._na33.cyhd.cychp.txt",
  "HV80R" = "/mnt/projects/p2ry8-crlf2/data/snparray/St.Anna_HV80-R.NA33.cyhd.cychp.txt"
)

regions = data.frame(patient=character(0), name=character(0), chr=character(0), start=numeric(0), end=numeric(0), genes=character(0), stringsAsFactors = F)

# PAR1 ---------------

regions = rbind(regions, setNames(data.frame("737D", "PAR1", "X", 1039250, 1984490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("737R", "PAR1", "X", 1039250, 1984490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("737R3", "PAR1", "X", 1039250, 1984490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))

#regions = rbind(regions, setNames(data.frame("715D", "PAR1", "X", 1039250, 1984490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))
#regions = rbind(regions, setNames(data.frame("715R", "PAR1", "X", 1039250, 1984490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("715R3", "PAR1", "X", 1039250, 1984490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("92D", "PAR1", "X", 1039250, 1984490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("92R", "PAR1", "X", 1039250, 1984490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("839D", "PAR1", "X", 1039250, 1984490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("839R", "PAR1", "X", 1039250, 1984490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("948D", "PAR1", "X", 1039250, 1984490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("HV80R", "PAR1", "X", 1039250, 1984490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))

# IKZF1 ---------------

regions = rbind(regions, setNames(data.frame("737D", "IKZF1", "7", 49925794, 50867611, "IKZF1", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("737R", "IKZF1", "7", 49925794, 50867611, "IKZF1", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("737R3", "IKZF1", "7", 49925794, 50867611, "IKZF1", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("715R3", "IKZF1", "7", 49925794, 50867611, "IKZF1", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("92D", "IKZF1", "7", 49925794, 50867611, "IKZF1", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("92R", "IKZF1", "7", 49925794, 50867611, "IKZF1", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("839D", "IKZF1", "7", 49925794, 50867611, "IKZF1", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("839R", "IKZF1", "7", 49925794, 50867611, "IKZF1", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("948D", "IKZF1", "7", 49925794, 50867611, "IKZF1", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("HV80R", "IKZF1", "7", 49925794, 50867611, "IKZF1", stringsAsFactors = F), names(regions)))

# PAX5 ---------------

regions = rbind(regions, setNames(data.frame("737D", "PAX5", "9", 36430065, 37406492, "PAX5", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("737R", "PAX5", "9", 36430065, 37406492, "PAX5", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("737R3", "PAX5", "9", 36430065, 37406492, "PAX5", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("715R3", "PAX5", "9", 36430065, 37406492, "PAX5", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("92D", "PAX5", "9", 36430065, 37406492, "PAX5", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("92R", "PAX5", "9", 36430065, 37406492, "PAX5", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("839D", "PAX5", "9", 36430065, 37406492, "PAX5", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("839R", "PAX5", "9", 36430065, 37406492, "PAX5", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("948D", "PAX5", "9", 36430065, 37406492, "PAX5", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("HV80R", "PAX5", "9", 36430065, 37406492, "PAX5", stringsAsFactors = F), names(regions)))

# CDKN2A ---------------

regions = rbind(regions, setNames(data.frame("737D", "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("737R", "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("737R3", "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("715R3", "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("92D", "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("92R", "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("839D", "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("839R", "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("948D", "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))

regions = rbind(regions, setNames(data.frame("HV80R", "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))


cols <- c("darkblue", "blue", "darkgray", "red", "darkred", "darkred", "darkred")
par(pch=19)

# collect data and build tracks; parallelize to speed up

tracks.all <- foreach(i=1:nrow(regions), .verbose = FALSE) %dopar% {
  tracks <- list()
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
  tracks[["ideogram"]] <- IdeogramTrack(genome = "hg19", chromosome = paste0("chr", region$chr), showId=F, showTitle=T, 
                                  name=paste0("chr", region$chr), fontcolor.title="black", background.title="white", rotation.title=0, fontsize=40)
  tracks[["axis"]] <- GenomeAxisTrack(labelPos="above", name=paste0("chr", region$chr), lwd=0.5, cex=0.25, col="black", fontcolor="black", fontface=2)
  tracks[["genes"]] <- BiomartGeneRegionTrack(genome = "hg19", chromosome = paste0("chr", region$chr), 
                                      start = region$start, end = region$end, name = "", showId=TRUE, 
                                      stacking="squish", collapseTranscripts = "longest",
                                      fontsize=7, fontcolor.group="black", background.title="white")
  tracks[["genes"]]@range <- tracks[["genes"]]@range[tracks[["genes"]]@range$symbol %in% unlist(strsplit(region$genes, ","))]
  tracks[["array"]] <- DataTrack(lrr.gr, name = paste(region$patient, "Array"), type = c("p", "g"), groups=c("Gain", "Loss", "Normal"), col=c("red", "blue", "darkgray"), 
                         cex=0.2, lty.grid=2, v=0, ylim=c(-1.8, 1.8), alpha.title = 1, alpha=0.9, fontsize=40, col.title="black", col.axis="black")
  tracks[["sep"]] <- DataTrack(ylim=c(1,1), showTitle=F, showAxis=F, background.panel = "lightgray")
  tracks[["wes"]] <- DataTrack(wes.gr, name = paste(region$patient, "WES"), type = c("p", "g"), groups=c("Gain", "Loss", "Normal"), col=c("darkgray", "blue", "red"), 
                         cex=0.4, lty.grid=2, v=0, ylim=c(0.1,1.9), alpha.title = 1, alpha=0.9, fontsize=40, col.title="black", col.axis="black")
  tracks
#  plot(
#    x=fit[[region$patient]][[region$chr]], 
#    main=paste(region$patient, region$name), 
#    cex=0.7, 
#    cex.main=1, 
#    cex.axis=1, 
#    ylim=c(0, 4), 
#    xlim=c(region$start, region$end), 
#    xlab=paste("Chromosome", region$chr), 
#    ylab="Copy number", 
#    col=cols,
#    show.legend=F,
#    panel.first=abline(h=1:4, lty=2, col="#dddddd", lwd=1)
#  )
}

# plot to PDF

pdf("/mnt/projects/p2ry8-crlf2/results/figures/exomeCopy.pdf", width=5.5/2.54, height=4/2.54)
for (i in 1:length(tracks.all)) {
  print(paste("Plotting", regions[i, "patient"], regions[i, "name"], "..."))
  plotTracks(
    trackList = tracks.all[[i]], 
    sizes     = c(0.09, 0.195, 0.09, 0.31, 0.005, 0.31),
    from      = regions[i, "start"], 
    to        = regions[i, "end"],
    cex.title = 0.1, cex.axis=0.08
  )
}
dev.off()
