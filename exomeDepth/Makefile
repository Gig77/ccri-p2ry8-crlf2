export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
.SECONDEXPANSION:

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PATIENTS_TEST = 92
PATIENTS_MATCHED = 108 92 737 839 B36 BB16 DL2 GI8 GI13 HV57 HV80 LU3 MA5 N7 S23 SN18 DS10898 VS14645 SE15285 BJ17183 KE17247 AL9890 GL11356
PATIENTS_DIA_ONLY = 242 360 365 379 400 506 769 802 833 841 887 903 948 957 961 1060 1066 1089 HW11537 KT14158 TL14516
PATIENTS_REL2 = 108 737 AL9890
PATIENTS_REL3 = 715 737 S23
PATIENTS_HD = 715 545 460 564 # sequenced with different protocol (TruSeq); CN for these samples is determined within the HD cohort
PATIENTS_XENO = m1963-545-rel m1964-545-rel m1957-715-rel m1977-G-dia m1967-Y-rel m1966-Y-dia m1035-108-dia m1060-108-rel m252-379-dia m1041-737-dia m1069-737-rel m247-833-dia m248-841-dia m1059-92-dia m1037-839-dia

all: allsamples.cnvs.exomeDepth.tsv

counts.RData: /mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_S23C_realigned.bam /mnt/projects/generic/scripts/exomeDepth/get-bam-counts.R
	Rscript /mnt/projects/generic/scripts/exomeDepth/get-bam-counts.R 2>&1 | $(LOG)

allsamples.cnvs.exomeDepth.tsv: $(foreach P, $(PATIENTS_MATCHED), $PC.cnvs.tsv $PD.cnvs.tsv $PR.cnvs.tsv) \
				   			    $(foreach P, $(PATIENTS_DIA_ONLY), $PC.cnvs.tsv $PD.cnvs.tsv) \
					 			$(foreach P, $(PATIENTS_REL2), $PR2.cnvs.tsv) \
					 			$(foreach P, $(PATIENTS_REL3), $PR3.cnvs.tsv) \
					 			$(foreach P, $(PATIENTS_XENO), $P.cnvs.tsv) \
					 			/mnt/projects/p2ry8-crlf2/scripts/exomeDepth/merge-cnvs.R
	Rscript /mnt/projects/p2ry8-crlf2/scripts/exomeDepth/merge-cnvs.R 2>&1 | $(LOG)
	mv $@.part $@

%.cnvs.tsv: counts.RData /mnt/projects/p2ry8-crlf2/scripts/exomeDepth/get-cnvs.R
	Rscript /mnt/projects/p2ry8-crlf2/scripts/exomeDepth/get-cnvs.R --sample $* 2>&1 | $(LOG)
	mv $*.cnvs.tsv.part $*.cnvs.tsv
	if [[ -a $*.exomeDepth.pdf.part ]]; then mv $*.exomeDepth.pdf.part $*.exomeDepth.pdf; fi;

#allpatients.filtered-segments.exomeCopy.tsv: allpatients.compiled-segments.exomeCopy.tsv /mnt/projects/generic/scripts/exomeCopy/filter-segments.R
#	Rscript /mnt/projects/generic/scripts/exomeCopy/filter-segments.R 2>&1 | $(LOG)
#	mv $@.part $@
	 
