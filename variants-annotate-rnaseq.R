#install.packages("doSNOW")
library(foreach)
library(doSNOW)
cl <- makeCluster(detectCores(), outfile="")
registerDoSNOW(cl)

s <- read.delim("/mnt/projects/ikaros/data/samples.csv", stringsAsFactors = F)
s <- s[s$Xeno != "yes",]
s$Patient[s$Patient=="AL"] <- "AL9890"

m <- read.delim("/mnt/projects/p2ry8-crlf2/results/filtered-variants.cosmic.tsv", stringsAsFactors = F)
m <- m[m$status != "REJECT",]
m <- m[m$var_type == "snp",]

totres <- foreach(i=1:nrow(s)) %dopar% {
  bam <- paste0("/mnt/projects/ikaros/results/anduril/execute/gsnap_", s$Alias[i], "/", s$Alias[i], ".gsnap.sorted.dupmarked.bam")
  sample <- ifelse(s$Timepoint[i]=="diagnosis", "rem_dia", "rem_rel")
  ms <- m[m$patient==s$Patient[i] & m$sample == sample,]
  
  res <- data.frame(patient=character(0), sample=character(0), chr=character(0), pos=integer(0), ref=character(0), alt=character(0), numtot=integer(0), numalt=integer(0))
  for (j in 1:nrow(ms)) {
    cmd <- paste0("/data_synology/software/samtools-0.1.19/samtools view -bh ", bam, " ", ms$chr[j], ":", ms$pos[j], "-", ms$pos[j]+1, " | /data_synology/software/samtools-0.1.19/samtools mpileup -q 1 -f /data_synology/anduril/docker-images/anduril-gsnap_2014_12_28-human_g1k_v37/db/human_g1k_v37/human_g1k_v37_etv6runx1.fasta - 2>/dev/null | grep ", ms$pos[j])
    r <- system(cmd, intern=T)
    if(length(r) > 0) {
      bases <- strsplit(r, "\t")[[1]][5]
      numalt <- sum(unlist(strsplit(bases, "")) %in% c(tolower(ms$alt[j]), toupper(ms$alt[j])))
      numtot <- numalt + sum(unlist(strsplit(bases, "")) %in% c(".", ","))
    } else {
      numalt <- 0
      numtot <- 0
    }
    res <- rbind(res, data.frame(patient=s$Patient[i], sample=sample, chr=ms$chr[j], pos=ms$pos[j], ref=ms$ref[j], alt=ms$alt[j], numtot=numtot, numalt=numalt))
    print(paste(i, "of", nrow(s), j, "of", nrow(ms), s$Alias[i], s$Patient[i], s$Timepoint[i], paste0(ms$chr[j], ":", ms$pos[j]), paste0(ms$ref[j], ">", ms$alt[j]), ms$dp_leu_ref[j], ms$dp_leu_var[j], ms$freq_leu[j], numtot, numalt))
  }
  res
}

# combine results into single dataframe
totres <- do.call("rbind", totres)
stopCluster(cl)

# merge with mutations
ma <- read.delim("/mnt/projects/p2ry8-crlf2/results/filtered-variants.cosmic.tsv", stringsAsFactors = F)
merged <- merge(ma, totres, by=c("patient", "sample", "chr", "pos", "ref", "alt"), all.x=T)
write.table(merged, "/mnt/projects/p2ry8-crlf2/results/filtered-variants.cosmic.rnaseq.tsv", col.names=T, row.names=F, sep="\t", quote=F, na="")
