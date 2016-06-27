library(foreach)
library(doMC)
library(Gviz)
library(GenomicRanges)
library(exomeCopy)

registerDoMC(40) # number of cores for parallel execution

rm(list=ls())

cychp = list(
  "108R2"    = "/mnt/projects/p2ry8-crlf2/data/snparray/108_R_4.cyhd.cychp.txt",
  "737D"     = "/mnt/projects/p2ry8-crlf2/data/snparray/A1311-15.CEL._2007-2323.na33.cyhd.cychp.txt",
  "737R"     = "/mnt/projects/p2ry8-crlf2/data/snparray/2009-2490_A1311-16.CEL._na33.cyhd.cychp.txt",
  "737R3"    = "/mnt/projects/p2ry8-crlf2/data/snparray/2012-2407._na33.cyhd.cychp.txt",
#  "715D"     = "/mnt/projects/p2ry8-crlf2/data/snparray/2007-0966_A1311-12.CEL.NA33.cyhd.cychp.txt",   # low-quality WES data
#  "715R"     = "/mnt/projects/p2ry8-crlf2/data/snparray/2008-2775_A1311-13.CEL..NA33.cyhd.cychp.txt",  # low-quality WES data
  "715R3"    = "/mnt/projects/p2ry8-crlf2/data/snparray/715_RR._na33.cyhd.ND.cychp.txt",
  "92D"      = "/mnt/projects/p2ry8-crlf2/data/snparray/94-1115_92_D._na33.cyhd.cychp.txt",
  "92R"      = "/mnt/projects/p2ry8-crlf2/data/snparray/98-0366_92_R._na33.cyhd.cychp.txt",
  "839D"     = "/mnt/projects/p2ry8-crlf2/data/snparray/2009-1595_839_D._na33.cyhd.cychp.txt",
  "839R"     = "/mnt/projects/p2ry8-crlf2/data/snparray/2013-0667_839_R._na33.cyhd.cychp.txt",
  "948D"     = "/mnt/projects/p2ry8-crlf2/data/snparray/2010-3831._na33.cyhd.cychp.txt",
  "HV80R"    = "/mnt/projects/p2ry8-crlf2/data/snparray/St.Anna_HV80-R.NA33.cyhd.cychp.txt",
  "242D"     = "/mnt/projects/p2ry8-crlf2/data/snparray/99-2714_A1311-08.CEL.na33.cyhd.cychp.txt",
  "379D"     = "/mnt/projects/p2ry8-crlf2/data/snparray/379D.CEL.cyhd.cychp.txt", # 2001-3936
  "957D"     = "/mnt/projects/p2ry8-crlf2/data/snparray/957D.CEL.cyhd.cychp.txt", # 2011-0318
  "961D"     = "/mnt/projects/p2ry8-crlf2/data/snparray/961D.CEL.cyhd.cychp.txt", # 2011-0419
  "KT14158D" = "/mnt/projects/p2ry8-crlf2/data/snparray/KT14158-Dg_St.Anna.NA33.cyhd.cychp.txt"
)

regions = data.frame(patient=character(0), name=character(0), chr=character(0), start=numeric(0), end=numeric(0), genes=character(0), stringsAsFactors = F)

# PAR1 ---------------

for (s in names(cychp)) {
  regions = rbind(regions, setNames(data.frame(s, "PAR1", "X", 1039250, 1984490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))
}

# IKZF1 ---------------

for (s in c(names(cychp), "N7D", "N7R", "DL2D", "DL2R", "HV57D", "HV57R")) {
  regions = rbind(regions, setNames(data.frame(s, "IKZF1", "7", 49925794, 50867611, "IKZF1", stringsAsFactors = F), names(regions)))
  regions = rbind(regions, setNames(data.frame(s, "IKZF1", "7", 50300000, 50500000, "IKZF1", stringsAsFactors = F), names(regions))) # zoom in
}

# PAX5 ---------------

for (s in names(cychp)) {
  regions = rbind(regions, setNames(data.frame(s, "PAX5", "9", 36430065, 37406492, "PAX5", stringsAsFactors = F), names(regions)))
}

# CDKN2A ---------------

for (s in names(cychp)) {
  regions = rbind(regions, setNames(data.frame(s, "CDKN2A", "9", 21519360, 22469828, "CDKN2A", stringsAsFactors = F), names(regions)))
}

# get array CNV segmented data
cnv.array <- read.delim("/mnt/projects/p2ry8-crlf2/results/cnvs.snp_arrays.txt", stringsAsFactor=F)
cnv.array$seqnames <- paste0("chr", cnv.array$seqnames)
names(cnv.array)[names(cnv.array)=="copy.count"] <- "Value"

cols <- c("darkblue", "blue", "darkgray", "red", "darkred", "darkred", "darkred")
par(pch=19)

# collect data and build tracks; parallelize to speed up

#tracks.all <- foreach(i=5, .verbose = TRUE) %dopar% {
tracks.all <- foreach(i=1:nrow(regions), .verbose = TRUE) %dopar% {
  tracks <- list()
  region <- regions[i,]

  # get array data
  if (!is.null(cychp[[region$patient]])) {
    cmd <- paste0("cat ", cychp[[region$patient]], " | perl /mnt/projects/p2ry8-crlf2/scripts/cychp2csv.pl ", region$chr, ":", region$start, "-", region$end)
    print(cmd)
    array <- read.table(pipe(cmd), header=T, stringsAsFactors = F)
    array$Chromosome <- paste0("chr", array$Chromosome)
    lrr.gr <- makeGRangesFromDataFrame(array[array$Type=="LRR",c("Chromosome", "Start", "Value")], start.field = "Start", end.field = "Start", keep.extra.columns = T)
  
    # prefer CN calls from externally provided file, otherwise use CN segments stored inside .cychp file
    if (region$patient %in% cnv.array$sample) {  
      cn.gr <- makeGRangesFromDataFrame(cnv.array[cnv.array$sample == region$patient & cnv.array$Value != 2, c("seqnames", "start", "end", "Value")], start.field = "start", end.field = "end", keep.extra.columns = T)
      
      # fill in diploid gaps in plotting region
      diff <- setdiff(GRanges(paste0("chr", region$chr), IRanges(region$start, region$end)), cn.gr)
      if (length(diff) > 0) {
        diff$Value <- 2
        cn.gr <- c(cn.gr, diff)
      }
    } else {
      cn.gr <- makeGRangesFromDataFrame(array[array$Type=="CN",c("Chromosome", "Start", "End", "Value")], start.field = "Start", end.field = "End", keep.extra.columns = T)
    }
    
    # encode CN state into LRR
    o <- findOverlaps(lrr.gr, cn.gr)
    states <- cn.gr[o@subjectHits]$Value
    lrr.gr <- GRanges(seqnames=seqnames(lrr.gr), ranges=ranges(lrr.gr), 
                      mcols=data.frame(
                        Gain   = as.numeric(ifelse(states >  2, lrr.gr$Value, NA)),
                        Loss   = as.numeric(ifelse(states <  2, lrr.gr$Value, NA)),
                        Normal = as.numeric(ifelse(states == 2, lrr.gr$Value, NA))
                      ))
  } else {
    lrr.gr <- GRanges(seqnames=region$chr, ranges=IRanges(1,2), mcols=data.frame(Gain=2, Loss=2, Normal=2))
  }

  # get WES data
  load(paste0("/mnt/projects/p2ry8-crlf2/results/exomeCopy/", region$patient, ".exomeCopy.fit.RData"))
  states <- fit[[1]][[region$chr]]@fx.par$S
  path <- as.vector(fit[[1]][[region$chr]]@path)
  ratio <- fit[[1]][[region$chr]]@O.norm
  wes.gr <- GRanges(seqnames=paste0("chr", region$chr), 
                    ranges=fit[[1]][[region$chr]]@ranges[[1]], 
                    mcols=data.frame(
                      Normal=as.numeric(ifelse(states[path] == 2, ratio, NA)),
                      Loss=as.numeric(ifelse(states[path] < 2, ratio, NA)),
                      Gain=as.numeric(ifelse(states[path] > 2, ratio, NA))
                    ))

  # configure tracks
  patient.label <- region$patient
  if (patient.label == "737R3") patient.label <- "737RR"
  if (patient.label == "108R2") patient.label <- "108RR"
  
  tracks[["ideogram"]] <- IdeogramTrack(genome = "hg19", chromosome = paste0("chr", region$chr), showId=F, showTitle=T, 
                                  name=paste0("chr", region$chr), fontcolor.title="black", background.title="white", rotation.title=0, fontsize=40)
  tracks[["axis"]] <- GenomeAxisTrack(labelPos="above", name=paste0("chr", region$chr), lwd=0.5, cex=0.25, col="black", fontcolor="black", fontface=2)
  tracks[["genes"]] <- BiomartGeneRegionTrack(genome = "hg19", chromosome = paste0("chr", region$chr), 
                                      start = region$start, end = region$end, name = "", showId=TRUE, 
                                      stacking="squish", collapseTranscripts = "longest",
                                      fontsize=7, fontcolor.group="black", background.title="white")
  tracks[["genes"]]@range <- tracks[["genes"]]@range[tracks[["genes"]]@range$symbol %in% unlist(strsplit(region$genes, ","))]
  tracks[["array"]] <- DataTrack(lrr.gr, name = paste(patient.label, "Array"), type = c("p", "g"), groups=c("Gain", "Loss", "Normal"), col=c("red", "blue", "darkgray"), 
                       cex=0.2, lty.grid=2, v=0, ylim=c(-1.8, 1.8), alpha.title = 1, alpha=0.9, fontsize=40, col.title="black", col.axis="black")
  tracks[["sep"]] <- DataTrack(ylim=c(1,1), showTitle=F, showAxis=F, background.panel = "lightgray")
  tracks[["wes"]] <- DataTrack(wes.gr, name = paste(patient.label, "WES"), type = c("p", "g"), groups=c("Gain", "Loss", "Normal"), col=c("darkgray", "blue", "red"), 
                         cex=0.4, lty.grid=2, v=0, ylim=c(0.1,1.9), alpha.title = 1, alpha=0.9, fontsize=40, col.title="black", col.axis="black")
  tracks
#  plot(
#    x=fit[[region$patient]][[region$chr]], 
#    main=paste(patient.label, region$name), 
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
