export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
.SECONDEXPANSION:

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PATIENTS_TEST = 92
PATIENTS_MATCHED = 108 839 92 B36 BB16 GI13 HV57 HV80 N7 S23 737
PATIENTS_DIA_ONLY = 242 360 365 379 400 506 769 833 948
PATIENTS_REL2 = 108 737 
PATIENTS_REL3 = 715 
PATIENTS_EXCLUDED = LU3 SN18

all: circos/allpatients.relapse.circos.png circos/allpatients.diagnosis.circos.png

counts.p2ry8.RData: ~/p2ry8-crlf2/data/bam/108C.duplicate_marked.realigned.recalibrated.bam
	Rscript ~/generic/scripts/exomeCopy/get-counts.R

counts.bg.p2ry8.RData: counts.p2ry8.RData ~/generic/scripts/exomeCopy/calc-background.R
	Rscript ~/generic/scripts/exomeCopy/calc-background.R

%.combined.pdf %.compiled-segments.tsv: counts.bg.p2ry8.RData ~/generic/scripts/exomeCopy/get-cnv.R
	Rscript ~/generic/scripts/exomeCopy/get-cnv.R --sample S$*

allpatients.compiled-segments.exomeCopy.tsv: $(foreach P, $(PATIENTS_MATCHED), $PC.compiled-segments.tsv $PD.compiled-segments.tsv $PR.compiled-segments.tsv) \
                                             $(foreach P, $(PATIENTS_DIA_ONLY), $PC.compiled-segments.tsv $PD.compiled-segments.tsv) \
                                             $(foreach P, $(PATIENTS_REL2), $PR2.compiled-segments.tsv) \
                                             ~/generic/scripts/exomeCopy/merge-cnvs.R
	Rscript ~/generic/scripts/exomeCopy/merge-cnvs.R

circos/%.somatic.circos.tsv: allpatients.compiled-segments.exomeCopy.tsv ~/p2ry8-crlf2/scripts/exomeCopy/get-circos-somatic.R
	Rscript ~/p2ry8-crlf2/scripts/exomeCopy/get-circos-somatic.R --segments-file $< --output-file $@.part --tumor S$* --normal wurscht
	mv $@.part $@ 

circos/allpatients.diagnosis.circos.png: $(foreach P, $(PATIENTS_MATCHED) $(PATIENTS_DIA_ONLY), circos/$PD.somatic.circos.tsv) ~/p2ry8-crlf2/scripts/exomeCopy/template.circos.conf ~/p2ry8-crlf2/scripts/exomeCopy/make-circos-conf.pl ~/p2ry8-crlf2/scripts/exomeCopy/gene-labels.txt
	rm -f circos/diagnosis.circos.conf
	perl ~/p2ry8-crlf2/scripts/exomeCopy/make-circos-conf.pl \
		--template ~/p2ry8-crlf2/scripts/exomeCopy/template.circos.conf \
		--data-dir ~/p2ry8-crlf2/results/exomeCopy/circos \
		--file-pattern D.somatic.circos.tsv \
		--gfx-output-file $@.part \
		> circos/diagnosis.circos.conf
	~/tools/circos-0.64/bin/circos --conf circos/diagnosis.circos.conf
	mv $@.part.png $@
	
circos/allpatients.relapse.circos.png: $(foreach P, $(PATIENTS_MATCHED), circos/$PR.somatic.circos.tsv) ~/p2ry8-crlf2/scripts/exomeCopy/template.circos.conf ~/p2ry8-crlf2/scripts/exomeCopy/make-circos-conf.pl ~/p2ry8-crlf2/scripts/exomeCopy/gene-labels.txt
	rm -f circos/relapse.circos.conf
	perl ~/p2ry8-crlf2/scripts/exomeCopy/make-circos-conf.pl \
		--template ~/p2ry8-crlf2/scripts/exomeCopy/template.circos.conf \
		--data-dir ~/p2ry8-crlf2/results/exomeCopy/circos \
		--file-pattern R.somatic.circos.tsv \
		--gfx-output-file $@.part \
		> circos/relapse.circos.conf
	~/tools/circos-0.64/bin/circos --conf circos/relapse.circos.conf
	mv $@.part.png $@