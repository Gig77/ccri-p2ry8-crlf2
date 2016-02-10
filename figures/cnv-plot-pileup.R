library(Gviz)
library(GenomicRanges)

# regions
regions = data.frame(name=character(0), chr=character(0), start=numeric(0), end=numeric(0), genes=character(0), stringsAsFactors = F)
regions = rbind(regions, setNames(data.frame("IKZF1", "chr7", 50308167, 50497707, "IKZF1", stringsAsFactors = F), names(regions)))
#regions = rbind(regions, setNames(data.frame("PAR1", "chrX", 539250, 2484490, "P2RY8,CRLF2", stringsAsFactors = F), names(regions)))

# samples
samples = data.frame(id=character(0), region=character(0), source=character(0), stringsAsFactors = F)
samples = rbind(samples, setNames(data.frame("N7D", "IKZF1", "PCR", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("N7R", "IKZF1", "WES", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("HV57D", "IKZF1", "MLPA", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("HV57R", "IKZF1", "MLPA", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("GI8R", "IKZF1", "MLPA", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("DL2D", "IKZF1", "MLPA", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("DL2R", "IKZF1", "MLPA", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("108D", "IKZF1", "MLPA", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("108R", "IKZF1", "MLPA", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("108R2", "IKZF1", "WES", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("DS10898D", "IKZF1", "MLPA", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("DS10898R", "IKZF1", "MLPA", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("KT14158D", "IKZF1", "Array", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("BJ17183D", "IKZF1", "WES", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("HV80D", "IKZF1", "MLPA", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("HV80R", "IKZF1", "Array", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("841D", "IKZF1", "WES", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("B36R", "IKZF1", "MLPA", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("92D", "IKZF1", "Array", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("92R", "IKZF1", "Array", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("GI13R", "IKZF1", "WES", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("737D", "IKZF1", "Array", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("737R", "IKZF1", "Array", stringsAsFactors = F), names(samples)))
samples = rbind(samples, setNames(data.frame("737R3", "IKZF1", "Array", stringsAsFactors = F), names(samples)))

#region <- "IKZF1"; chr <- "chr7"; start <- 50000000; end <- 50800000
#region <- "CDKN2AB"; chr <- "chr9"; start <- 21600000; end <- 22400000
#region <- "PAR1"; chr <- "chrX"; start <- 1000000; end <- 2000000
#region <- "PAX5"; chr <- "chr9"; start <- 36495626; end <- 37377382
#region <- "GABRB3"; chr <- "chr15"; start <- 25870574; end <- 27936343
#region <- "ACSM2A"; chr <- "chr16"; start <- 20318327; end <- 20643523
#region <- "BTG1"; chr <- "chr12"; start <- 91278699; end <- 92889529

# get WES CNVs
cnv.wes <- read.delim("/mnt/projects/p2ry8-crlf2/results/exomeCopy/allpatients.filtered-segments.exomeCopy.tsv", stringsAsFactors = F)
cnv.wes$seqnames <- paste0("chr", cnv.wes$seqnames)
cnv.wes$start <- cnv.wes$start-2000 ; cnv.wes$end <- cnv.wes$end+2000  # better align segments with exons of gene model
cnv.wes.gr <- makeGRangesFromDataFrame(cnv.wes, keep.extra.columns = T)

# get array CNVs
cnv.arr <- read.delim("/mnt/projects/p2ry8-crlf2/results/cnvs.snp_arrays.txt", stringsAsFactors = F)
cnv.arr$seqnames <- paste0("chr", cnv.arr$seqnames)
cnv.arr.gr <- makeGRangesFromDataFrame(cnv.arr, keep.extra.columns = T)

# MLPA CNVs
cnv.mlpa <- data.frame(sample.name=character(0), seqnames=character(0), start=numeric(0), end=numeric(0), copy.count=integer(0))
cnv.mlpa <- rbind(cnv.mlpa, setNames(data.frame("DL2D", "chr7", 50442236, 50475384, 1, stringsAsFactors = F), names(cnv.mlpa)))  # exon 4-8
cnv.mlpa <- rbind(cnv.mlpa, setNames(data.frame("DL2R", "chr7", 50442236, 50475384, 1, stringsAsFactors = F), names(cnv.mlpa)))  # exon 4-8
cnv.mlpa <- rbind(cnv.mlpa, setNames(data.frame("108D", "chr7", 50442236, 50475384, 1, stringsAsFactors = F), names(cnv.mlpa)))  # exon 4-8
cnv.mlpa <- rbind(cnv.mlpa, setNames(data.frame("108R", "chr7", 50442236, 50475384, 1, stringsAsFactors = F), names(cnv.mlpa)))  # exon 4-8
cnv.mlpa <- rbind(cnv.mlpa, setNames(data.frame("HV57D", "chr7", 50442236, 50462204, 1, stringsAsFactors = F), names(cnv.mlpa))) # exon 4-7
cnv.mlpa <- rbind(cnv.mlpa, setNames(data.frame("HV57R", "chr7", 50442236, 50462204, 1, stringsAsFactors = F), names(cnv.mlpa))) # exon 4-7
cnv.mlpa <- rbind(cnv.mlpa, setNames(data.frame("HV80D", "chr7", 50442236, 50607726, 1, stringsAsFactors = F), names(cnv.mlpa))) # exon 4-8 (incl. neighboring gene DDC)
cnv.mlpa <- rbind(cnv.mlpa, setNames(data.frame("HV80R", "chr7", 50442236, 50607726, 1, stringsAsFactors = F), names(cnv.mlpa))) # exon 4-8 (incl. neighboring gene DDC)
cnv.mlpa <- rbind(cnv.mlpa, setNames(data.frame("DS10898D", "chr7", 50342375, 50462204, 1, stringsAsFactors = F), names(cnv.mlpa))) # exon 1-7
cnv.mlpa <- rbind(cnv.mlpa, setNames(data.frame("DS10898R", "chr7", 50342375, 50462204, 1, stringsAsFactors = F), names(cnv.mlpa))) # exon 1-7
cnv.mlpa <- rbind(cnv.mlpa, setNames(data.frame("GI8R", "chr7", 50442236, 50462204, 1, stringsAsFactors = F), names(cnv.mlpa)))  # exon 4-7
cnv.mlpa <- rbind(cnv.mlpa, setNames(data.frame("B36R", "chr7", 50000000, 50475384, 1, stringsAsFactors = F), names(cnv.mlpa)))  # exon 1-8
cnv.mlpa <- rbind(cnv.mlpa, setNames(data.frame("KT14158D", "chr7", 50442236, 50475384, 1, stringsAsFactors = F), names(cnv.mlpa))) # exon 4-8
cnv.mlpa.gr <- makeGRangesFromDataFrame(cnv.mlpa, keep.extra.columns = T)

# PCR CNVs
cnv.pcr <- data.frame(sample.name=character(0), seqnames=character(0), start=numeric(0), end=numeric(0), copy.count=integer(0))
cnv.pcr <- rbind(cnv.pcr, setNames(data.frame("N7D", "chr7", 50442236, 50462204, 1, stringsAsFactors = F), names(cnv.pcr))) # exon 4-7
cnv.pcr.gr <- makeGRangesFromDataFrame(cnv.pcr, keep.extra.columns = T)

#d <- rbind(d,data.frame(source="MLPA", seqnames="7", start=50435703, end=50459561, width=23859, strand="*", sample.name="108R", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="IKZF1", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
#d <- rbind(d,data.frame(source="MLPA", seqnames="7", start=50435703, end=50459561, width=23859, strand="*", sample.name="DL2R", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="IKZF1", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
#d <- rbind(d,data.frame(source="MLPA", seqnames="7", start=50358657, end=50459561, width=100905, strand="*", sample.name="DS10898D", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="IKZF1", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 
#d <- rbind(d,data.frame(source="MLPA", seqnames="7", start=50358657, end=50459561, width=100905, strand="*", sample.name="GI8R", copy.count=1, log.odds=NA, nranges=NA, targeted.bp=NA, genes="IKZF1", overlap.samples=NA, overlap.count=NA, overlap.count.tumor=NA)) 

for (i in 1:nrow(regions)) {
  region <- regions[i,]
  samples.region <- samples[samples$region==region$name,]
  
  ideogram.track <- IdeogramTrack(genome = "hg19", chromosome = region$chr, showId=F, showTitle=T, 
                                  name=region$chr, fontcolor.title="black", background.title="white", rotation.title=0, fontsize=25)
  axis.track <- GenomeAxisTrack(labelPos="above", name=region$chr, lwd=0.5, cex=0.25, col="black", fontcolor="black", fontface=2)
  gene.track <- BiomartGeneRegionTrack(genome = "hg19", chromosome = region$chr, start = region$start, end = region$end, 
                                      name = "", showId=TRUE, stacking="squish", collapseTranscripts = "longest",
                                      fontsize=7, fontcolor.group="black", background.title="white")
  gene.track@range <- gene.track@range[gene.track@range$symbol %in% unlist(strsplit(region$genes, ","))]

  cntracks <- list()
  for (j in 1:nrow(samples.region)) {
    sample <- samples.region[j,]
    print(paste0("Adding track for sample ", sample$id, " (", sample$source, ")"))
    if (sample$source == "Array") {
      grs <- cnv.arr.gr[seqnames(cnv.arr.gr)==region$chr & cnv.arr.gr$sample==sample$id & cnv.arr.gr$type %in% c("Loss", "Gain"),]
    } else if (sample$source == "PCR") {
      grs <- cnv.pcr.gr[seqnames(cnv.pcr.gr)==region$chr & cnv.pcr.gr$sample.name==sample$id,]
    } else if (sample$source == "MLPA") {
      grs <- cnv.mlpa.gr[seqnames(cnv.mlpa.gr)==region$chr & cnv.mlpa.gr$sample.name==sample$id,]
    } else {
      grs <- cnv.wes.gr[seqnames(cnv.wes.gr)==region$chr & cnv.wes.gr$sample.name==sample$id,]
    }
    detected.by <- character(0)
    #if (length(findOverlaps(grs[grs$copy.count<2], cnv.wes.gr[cnv.wes.gr$sample.name==sample$id & cnv.wes.gr$copy.count<2])) > 0) detected.by <- c(detected.by, "E")
    #if (length(findOverlaps(grs[grs$copy.count<2], cnv.arr.gr[cnv.arr.gr$sample==sample$id & cnv.arr.gr$copy.count<2])) > 0) detected.by <- c(detected.by, "A")
    #if (length(findOverlaps(grs[grs$copy.count<2], cnv.pcr.gr[cnv.pcr.gr$sample.name==sample$id & cnv.pcr.gr$copy.count<2])) > 0) detected.by <- c(detected.by, "P")
    #if (length(findOverlaps(grs[grs$copy.count<2], cnv.mlpa.gr[cnv.mlpa.gr$sample.name==sample$id & cnv.mlpa.gr$copy.count<2])) > 0) detected.by <- c(detected.by, "M")
    cn.track <- AnnotationTrack(grs, name=paste(sample$id, paste(detected.by, collapse="/")), stackHeight=1, min.height=1, rot.title=0, 
                                col="white", lwd.border=0, fontsize=25, col.title="black", col.axis="black")
    if (length(grs) > 0) feature(cn.track) <- paste0("CN", grs$copy.count) 
    cntracks <- c(cntracks, cn.track)
  }

  pdf(paste0("/mnt/projects/p2ry8-crlf2/results/figures/cnv-allsamples.", region$name, ".pdf"), width=5.5/2.54, height=8/2.54)
  plotTracks(c(list(ideogram.track, axis.track, gene.track), cntracks), 
             sizes=c(0.045, 0.1, 0.05, rep((1-(0.045+0.1+0.05))/nrow(samples.region), nrow(samples.region))),
             from=region$start, to=region$end, title.width=1.8, 
             cex.title = 0.15, cex.axis=0.1,
             CN0="darkblue", "CN1"="blue", "CN2"="white", "CN3"="red", "CN4"="darkred")
  dev.off()
}
