export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
.SECONDEXPANSION:

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PATIENTS_TEST = 92
PATIENTS_MATCHED = 108 92 737 839 B36 BB16 DL2 GI8 GI13 HV57 HV80 LU3 MA5 N7 S23 SN18 DS10898 VS14645 SE15285 BJ17183 KE17247 AL9890 GL11356
PATIENTS_DIA_ONLY = 242 360 365 379 400 506 769 802 833 841 887 903 948 957 961 1060 1066 1089 HW11537 KT14158 TL14516
PATIENTS_REL2 = 108 737
PATIENTS_REL3 = 715 737
PATIENTS_HD = 715 545 460 564 # sequenced with different protocol (TruSeq); CN for these samples is determined within the HD cohort
PATIENTS_XENO = m1963-545-rel m1964-545-rel m1957-715-rel m1977-G-dia m1967-Y-rel

all: circos/allpatients.relapse.circos.png circos/allpatients.diagnosis.circos.png allpatients.compiled-segments.exomeCopy.tsv allpatients.filtered-segments.exomeCopy.tsv $(foreach P, $(PATIENTS_XENO), $P.compiled-segments.tsv) \

counts.RData: /mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_S23C_realigned.bam /mnt/projects/generic/scripts/exomeCopy/get-counts.R
	Rscript /mnt/projects/generic/scripts/exomeCopy/get-counts.R

counts.bg.RData: counts.RData /mnt/projects/generic/scripts/exomeCopy/calc-background.R counts.RData
	Rscript /mnt/projects/generic/scripts/exomeCopy/calc-background.R

%.combined.pdf %.compiled-segments.tsv: counts.bg.RData /mnt/projects/generic/scripts/exomeCopy/get-cnv.R
	Rscript /mnt/projects/generic/scripts/exomeCopy/get-cnv.R --sample $*

allpatients.compiled-segments.exomeCopy.tsv: $(foreach P, $(PATIENTS_MATCHED), $PC.compiled-segments.tsv $PD.compiled-segments.tsv $PR.compiled-segments.tsv) \
											 $(foreach P, $(PATIENTS_DIA_ONLY), $PC.compiled-segments.tsv $PD.compiled-segments.tsv) \
											 $(foreach P, $(PATIENTS_REL2), $PR2.compiled-segments.tsv) \
											 $(foreach P, $(PATIENTS_REL3), $PR3.compiled-segments.tsv) \
											 /mnt/projects/generic/scripts/exomeCopy/merge-cnvs.R
	Rscript /mnt/projects/generic/scripts/exomeCopy/merge-cnvs.R
	mv $@.part $@

allpatients.filtered-segments.exomeCopy.tsv: allpatients.compiled-segments.exomeCopy.tsv /mnt/projects/generic/scripts/exomeCopy/filter-segments.R
	Rscript /mnt/projects/generic/scripts/exomeCopy/filter-segments.R
	mv $@.part $@
	 
.PHONY: circos
circos: circos/allpatients.diagnosis.circos.png circos/allpatients.relapse.circos.png
 
circos/%.somatic.circos.tsv: allpatients.filtered-segments.exomeCopy.tsv /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/get-circos-somatic.R
	mkdir -p circos
	Rscript /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/get-circos-somatic.R --segments-file $< --output-file $@.part --tumor $* --normal wurscht
	mv $@.part $@ 

circos/allpatients.diagnosis.circos.png: $(foreach P, $(PATIENTS_MATCHED) $(PATIENTS_DIA_ONLY), circos/$PD.somatic.circos.tsv) /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/template.circos.conf /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/make-circos-conf.pl /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/gene-labels.txt /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/sample-order.dia.txt
	rm -f circos/diagnosis.circos.conf
	perl /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/make-circos-conf.pl \
		--template /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/template.circos.conf \
		--data-dir /mnt/projects/p2ry8-crlf2/results/exomeCopy/circos \
		--order-file /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/sample-order.dia.txt \
		--file-pattern D.somatic.circos.tsv \
		--gfx-output-file $@.part \
		> circos/diagnosis.circos.conf
	~/tools/circos-0.64/bin/circos --conf circos/diagnosis.circos.conf
	mv $@.part.png $@
	
circos/allpatients.relapse.circos.png: $(foreach P, $(PATIENTS_MATCHED), circos/$PR.somatic.circos.tsv) /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/template.circos.conf /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/make-circos-conf.pl /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/gene-labels.txt
	rm -f circos/relapse.circos.conf
	perl /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/make-circos-conf.pl \
		--template /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/template.circos.conf \
		--data-dir /mnt/projects/p2ry8-crlf2/results/exomeCopy/circos \
		--file-pattern R.somatic.circos.tsv \
		--gfx-output-file $@.part \
		> circos/relapse.circos.conf
	~/tools/circos-0.64/bin/circos --conf circos/relapse.circos.conf
	mv $@.part.png $@