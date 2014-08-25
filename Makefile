export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
.SECONDEXPANSION:

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PATIENTS_TEST = 92
PATIENTS_MATCHED = 108 839 92 B36 BB16 GI13 HV57 HV80 LU3 N7 S23 SN18 737
PATIENTS_DIA_ONLY = 242 360 365 379 400 506 769 833 948
PATIENTS_REL2 = 108 737 
PATIENTS_REL3 = 715 

#all: filtered-variants.cosmic.tsv filtered-variants.cosmic.normaf.tsv snpeff/mutect_737_rem_rel1.dbsnp.snpeff.dbNSFP.vcf snpeff/mutect_737_rem_rel2.dbsnp.snpeff.dbNSFP.vcf snpeff/mutect_737_rem_dia.dbsnp.snpeff.dbNSFP.vcf snpeff/indels_737_rem_rel1.indel.dbsnp.snpeff.dbNSFP.vcf snpeff/indels_737_rem_rel2.indel.dbsnp.snpeff.dbNSFP.vcf snpeff/indels_737_rem_dia.indel.dbsnp.snpeff.dbNSFP.vcf
all: filtered-variants.cosmic.normaf.tsv filtered-variants.cosmic.merged.tsv segmented_coverage/allpatients.segmented-coverage.pdf coverage/coverage-plots-exome.png picard

#-----------
# DOWNLOAD
#-----------
#wget -r -np -e robots=off --reject *.wig.txt,*_stats.out http://www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_zahyee0Zooxiphoogo5ooMee3mai8uy8/mutect_b37/ 
#wget -r -np -e robots=off --reject *.wig.txt,*_stats.out http://www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_zahyee0Zooxiphoogo5ooMee3mai8uy8/indels_b37/

~/p2ry8-crlf2/data/mutect_somatic_mutations_hg19/%_calls.vcf: ~/p2ry8-crlf2/data/mutect_somatic_mutations/%_calls.vcf
	cat $< | perl -ne 's/^([\dXY])/chr$$1/; s/^MT/chrM/; print $$_;' > $@.part
	mv $@.part $@

~/p2ry8-crlf2/data/mutect_somatic_indels_hg19/%_indel.vcf: ~/p2ry8-crlf2/data/mutect_somatic_indels/%_indel.vcf
	cat $< | perl -ne 's/^([\dXYM])/chr$$1/; print $$_;' > $@.part
	mv $@.part $@

#-----------
# REMISSION VARIANTS
#-----------
remission-variants.tsv: $(foreach P, $(PATIENTS_MATCHED), remission-variants/$P.remission-variants.snp.tsv remission-variants/$P.remission-variants.indel.tsv) \
					    $(foreach P, $(PATIENTS_DIA_ONLY), remission-variants/$P.remission-variants.snp.tsv remission-variants/$P.remission-variants.indel.tsv) \
					    ~/hdall/results/remission-variants.tsv
	cat $^ | ~/tools/lh3-sort/sort -k 2,2N -k 3,3n > remission-variants.tsv

remission-variants.tsv.gz: remission-variants.tsv
	bgzip -c $^ > $@.part
	mv $@.part $@

remission-variants.tsv.gz.tbi: remission-variants.tsv.gz
	~/tools/tabix-0.2.6/tabix $^ -s 2 -b 3 -e 3

remission-variants/%.remission-variants.snp.tsv: ~/p2ry8-crlf2/data/mutect_somatic_mutations_hg19/mutect_%_rem_dia_calls.vcf ~/p2ry8-crlf2/scripts/create-normal-vcf.pl
	cat $< | perl ~/p2ry8-crlf2/scripts/create-normal-vcf.pl \
		2>&1 1> $@.part | $(LOG)
	mv $@.part $@

remission-variants/%.remission-variants.indel.tsv: ~/p2ry8-crlf2/data/mutect_somatic_indels_hg19/indels_%_rem_dia_indel.vcf ~/p2ry8-crlf2/scripts/create-normal-vcf.pl
	cat $< | perl ~/p2ry8-crlf2/scripts/create-normal-vcf.pl \
		2>&1 1> $@.part | $(LOG)
	mv $@.part $@

#-----------	
# SNPEFF
#-----------	
snpeff/%.dbsnp.vcf: ~/p2ry8-crlf2/data/mutect_somatic_mutations_hg19/%_calls.vcf ~/generic/data/ncbi/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/generic/data/ncbi/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat $< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@

snpeff/%.indel.dbsnp.vcf: ~/p2ry8-crlf2/data/mutect_somatic_indels_hg19/%_indel.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/generic/data/ncbi/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat $< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@

snpeff/%.dbsnp.snpeff.vcf: snpeff/%.dbsnp.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -Xmx2g -jar snpEff.jar -v -lof hg19 -stats $(PWD)/snpeff/$*.snpeff.summary.html $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

snpeff/%.dbsnp.snpeff.dbNSFP.vcf: snpeff/%.dbsnp.snpeff.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar dbnsfp -v ~/generic/data/dbNSFP-2.1/dbNSFP2.1.txt $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@
	
#-------------
# CNV
#-------------

GENES_PAR1=PPP2R3B,SHOX,CRLF2,IL3RA,P2RY8,ASMT,DHR3X,ZBED1,CD99
GENES_IKZF1=VWC2,ZPBP,C7orf72,IKZF1,DDC,GRB10,COBL

.PHONY: coverage
coverage: coverage/allpatients.PAR1-X-60001-2699520.pdf \
		  coverage/allpatients.IKZF1-7-49578046-51601231.pdf \
		  coverage/allpatients.IKZF1_highres-7-50332565-50494236.pdf \
		  coverage/allpatients.IKZF2-2-211852462-217849831.pdf \
		  coverage/allpatients.IKZF3-17-37789665-38184756.pdf \
		  coverage/allpatients.PAX5-9-35800105-38084658.pdf \
		  coverage/allpatients.EBF1-5-155175854-161576670.pdf \
		  coverage/allpatients.ETV6-12-10309186-13646908.pdf \
		  coverage/allpatients.RUNX1-21-34971197-38263248.pdf \
		  coverage/allpatients.VPREB1-22-21866719-24317257.pdf \
		  coverage/allpatients.ERG-21-38349519-41562975.pdf \
		  coverage/allpatients.TP53-17-7440139-7721205.pdf \
		  coverage/allpatients.RB1-13-46907242-51051394.pdf \
		  coverage/allpatients.CDKN2AandB-9-20309360-23721195.pdf \
		  coverage/allpatients.CREBBP-16-3067833-4618486.pdf \
		  coverage/allpatients.MLL2-12-49295653-49584389.pdf \
		  coverage/allpatients.EZH2-7-148002918-149209618.pdf \
		  coverage/allpatients.NCOR1-17-15245025-16768909.pdf \
		  coverage/allpatients.TUSC3-8-12275019-19439531.pdf \
		  coverage/allpatients.WHSC1-4-1615018-2430195.pdf \
		  coverage/allpatients.NT5C2-10-104141780-105714035.pdf \
		  coverage/allpatients.LEF1-4-108313628-110095704.pdf \
		  coverage/allpatients.TCF3-19-1261454-2018719.pdf \
		  coverage/allpatients.BLNK-10-97442705-98708695.pdf \
		  coverage/allpatients.FOXO3A-6-108206749-109779719.pdf \
		  coverage/allpatients.FBXW7-4-151992600-154510682.pdf \
		  coverage/allpatients.CREG1-1-167207942-167849758.pdf \
		  coverage/allpatients.FLT3_PAN3-13-27955721-29384162.pdf \
		  coverage/allpatients.PDGFRB-5-149383560-149648693.pdf \
		  coverage/allpatients.STRN3-14-30942316-31995615.pdf \
		  coverage/allpatients.RANBP2-2-109007910-109671216.pdf \
		  coverage/allpatients.EPOR-19-11393474-11619196.pdf \
		  coverage/allpatients.SH2B3-12-111621200-112196224.pdf \
		  coverage/allpatients.HIST1H2BD_HIST1H1E-6-26033646-26316674.pdf \
		  coverage/allpatients.SPRED1-15-35957458-40803200.pdf \
		  coverage/allpatients.ADD3-10-111062840-113085360.pdf \
		  coverage/allpatients.ATP10A-15-24452756-27861838.pdf \
		  coverage/allpatients.NUP214-9-133609878-134537138.pdf
		  

segmented_coverage/allpatients.segmented-coverage.pdf: $(foreach P, $(PATIENTS_MATCHED), segmented_coverage/$PD.segmented-coverage.pdf segmented_coverage/$PR.segmented-coverage.pdf) \
								    				   $(foreach P, $(PATIENTS_DIA_ONLY), segmented_coverage/$PD.segmented-coverage.pdf)
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@
	rm $^

segmented_coverage/%D.segmented-coverage.pdf: segmented_coverage/%D.segmented-coverage.tsv segmented_coverage/%C.segmented-coverage.tsv ~/hdall/scripts/cnv/plot-segment-coverage.R
	Rscript ~/hdall/scripts/cnv/plot-segment-coverage.R \
		--patient $*D \
		--tumor $(word 1,$^) \
		--normal $(word 2,$^) \
		--circos circos/$*D.cnv.circos.tsv \
		--output $@.part
	mv $@.part $@

segmented_coverage/%R.segmented-coverage.pdf: segmented_coverage/%R.segmented-coverage.tsv segmented_coverage/%C.segmented-coverage.tsv ~/hdall/scripts/cnv/plot-segment-coverage.R
	Rscript ~/hdall/scripts/cnv/plot-segment-coverage.R \
		--patient $*R \
		--tumor $(word 1,$^) \
		--normal $(word 2,$^) \
		--circos circos/$*R.cnv.circos.tsv \
		--output $@.part
	mv $@.part $@

#segmented_coverage/%.segmented-coverage.tsv: ~/p2ry8-crlf2/data/bam/%.duplicate_marked.realigned.recalibrated.bam ~/hdall/scripts/cnv/get-segment-coverage.pl
segmented_coverage/%.segmented-coverage.tsv: ~/p2ry8-crlf2/data/bam/%.duplicate_marked.realigned.recalibrated.bam
	~/tools/samtools-0.1.19/samtools depth -Q 1 $< \
		| perl ~/hdall/scripts/cnv/get-segment-coverage.pl --sample $* --bin-size 250000 \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@
	
coverage/%.exon-coverage.tsv: ~/p2ry8-crlf2/data/bam/%.duplicate_marked.realigned.recalibrated.bam /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed ~/git/hdall/cnv/get-exon-coverage.pl
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-b /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed \
		$< \
	| perl ~/git/hdall/cnv/get-exon-coverage.pl \
		--exon-bed /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@

coverage/patient%.pdf: coverage/$$(word 1, $$(subst ., , %))D.exon-coverage.tsv coverage/$$(word 1, $$(subst ., , %))R.exon-coverage.tsv coverage/$$(word 1, $$(subst ., , %))C.exon-coverage.tsv ~/p2ry8-crlf2/scripts/cov-plot-region.R
	Rscript ~/p2ry8-crlf2/scripts/cov-plot-region.R \
		--patient $(word 1, $(subst ., , $*)) \
		--diagnosis $(word 1,$^) \
		--relapse $(word 2,$^) \
		--remission $(word 3,$^) \
		--output $@.part \
		--region-name $(word 1, $(subst -, , $(word 2, $(subst ., , $*)))) \
		--display-chrom $(word 2, $(subst -, , $*)) \
		--display-start $(word 3, $(subst -, , $*)) \
		--display-end $(word 4, $(subst -, , $*)) \
		$(if $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),--display-genes $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),)
	mv $@.part $@

coverage/patient%.diaonly.pdf: coverage/$$(word 1, $$(subst ., , %))D.exon-coverage.tsv coverage/$$(word 1, $$(subst ., , %))C.exon-coverage.tsv ~/p2ry8-crlf2/scripts/cov-plot-region.R
	Rscript ~/p2ry8-crlf2/scripts/cov-plot-region.R \
		--patient $(word 1, $(subst ., , $*)) \
		--diagnosis $(word 1,$^) \
		--remission $(word 2,$^) \
		--output $@.part \
		--region-name $(word 1, $(subst -, , $(word 2, $(subst ., , $*)))) \
		--display-chrom $(word 2, $(subst -, , $*)) \
		--display-start $(word 3, $(subst -, , $*)) \
		--display-end $(word 4, $(subst -, , $*)) \
		$(if $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),--display-genes $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),)
	mv $@.part $@

coverage/allpatients.%.pdf: $(foreach P, $(PATIENTS_MATCHED), coverage/patient$P.%.pdf) $(foreach P, $(PATIENTS_DIA_ONLY), coverage/patient$P.%.diaonly.pdf)
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@
	rm $^

coverage/coverage-plots-exome.png: $(foreach P, $(PATIENTS_MATCHED), coverage/$PD.coverage.bedtools.txt coverage/$PR.coverage.bedtools.txt coverage/$PC.coverage.bedtools.txt) \
               					   $(foreach P, $(PATIENTS_DIA_ONLY), coverage/$PD.coverage.bedtools.txt coverage/$PC.coverage.bedtools.txt) \
               					   ~/p2ry8-crlf2/scripts/plot-coverage.R
	Rscript ~/p2ry8-crlf2/scripts/plot-coverage.R

coverage/%.coverage.bedtools.txt: ~/p2ry8-crlf2/data/bam/%.duplicate_marked.realigned.recalibrated.bam ~/generic/data/illumina/nexterarapidcapture_exome_targetedregions.nochr.bed
	samtools view -bq 1 -F 0x400 $< | bedtools coverage -hist -abam - -b $(word 2, $^) | grep ^all > $@.part
	mv $@.part $@

segmented_coverage/%C.segmented-coverage.chr21.pdf: segmented_coverage/%C.segmented-coverage.tsv segmented_coverage/B36C.segmented-coverage.tsv ~/hdall/scripts/cnv/plot-segment-coverage-chr21.R
	Rscript ~/hdall/scripts/cnv/plot-segment-coverage-chr21.R \
		--patient $*C \
		--case $(word 1,$^) \
		--control segmented_coverage/B36C.segmented-coverage.tsv \
		--output $@.part
	mv $@.part $@

segmented_coverage/allpatients.segmented-coverage.chr21.pdf: $(foreach P, $(PATIENTS_MATCHED) $(PATIENTS_DIA_ONLY), segmented_coverage/$PC.segmented-coverage.chr21.pdf)
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@
	rm $^

#-------------
# final lists
#-------------
filtered-variants.tsv: $(foreach P, $(PATIENTS_MATCHED), filtered_variants/$P_rem_dia.snp.filtered.tsv filtered_variants/$P_rem_dia.indel.filtered.tsv filtered_variants/$P_rem_rel.snp.filtered.tsv filtered_variants/$P_rem_rel.indel.filtered.tsv) \
					   $(foreach P, $(PATIENTS_DIA_ONLY), filtered_variants/$P_rem_dia.snp.filtered.tsv filtered_variants/$P_rem_dia.indel.filtered.tsv) \
					   $(foreach P, $(PATIENTS_REL2), filtered_variants/$P_rem_rel2.snp.filtered.tsv filtered_variants/$P_rem_rel2.indel.filtered.tsv) \
#					   $(foreach P, $(PATIENTS_REL3), filtered_variants/$P_rem_rel3.snp.filtered.tsv filtered_variants/$P_rem_rel3.indel.filtered.tsv) \
					   ~/p2ry8-crlf2/scripts/filter-variants.pl 
	perl  ~/p2ry8-crlf2/scripts/filter-variants.pl --header 2>&1 1>$@.part | $(LOG)
	cat filtered_variants/*.filtered.tsv >> $@.part
	mv $@.part $@

filtered-variants.cosmic.tsv: filtered-variants.tsv ~/generic/data/cosmic/v67/CosmicMutantExport_v67_241013.tsv ~/hdall/scripts/annotate-cosmic.pl
	cat $(word 1,$^) | perl ~/hdall/scripts/annotate-cosmic.pl \
		--cosmic-mutation-file $(word 2,$^) \
		--only-confirmed \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@ 
	
filtered-variants.cosmic.normaf.tsv: filtered-variants.cosmic.tsv ~/hdall/results/cnv/hdall.cnv.tsv ~/hdall/scripts/normalize-af.pl
	cat filtered-variants.cosmic.tsv | perl ~/hdall/scripts/normalize-af.pl \
		--cnv-file ~/hdall/results/cnv/hdall.cnv.tsv \
		2>&1 1>$@.part | tee -a make.log
	mv $@.part $@ 

filtered-variants.cosmic.merged.tsv: filtered-variants.cosmic.tsv
	Rscript ~/p2ry8-crlf2/scripts/merge-filtered-variants.R 2>&1 | $(LOG)

filtered_variants/%.snp.filtered.tsv: snpeff/mutect_%.dbsnp.snpeff.dbNSFP.vcf ~/hdall/results/curated-recected-variants.tsv remission-variants.tsv.gz.tbi ~/p2ry8-crlf2/scripts/filter-variants.pl
	mkdir -p filtered_variants
	perl ~/p2ry8-crlf2/scripts/filter-variants.pl \
		--sample $* \
		--vcf-in $< \
		--variant-type snp \
		--vcf-out filtered_variants/$*.dbsnp.snpeff.dbNSFP.filtered.vcf \
		--rmsk-file ~/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file ~/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file ~/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file ~/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible ~/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro ~/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--rejected-variants-file ~/hdall/results/curated-recected-variants.tsv \
		--remission-variants-file remission-variants.tsv.gz \
		--evs-file ~/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		--min-num-rem-to-exclude 3 \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@

filtered_variants/%.indel.filtered.tsv: snpeff/indels_%.indel.dbsnp.snpeff.dbNSFP.vcf ~/hdall/results/curated-recected-variants.tsv remission-variants.tsv.gz.tbi ~/p2ry8-crlf2/scripts/filter-variants.pl
	mkdir -p filtered_variants
	perl ~/p2ry8-crlf2/scripts/filter-variants.pl \
		--sample $* \
		--vcf-in $< \
		--variant-type indel \
		--vcf-out filtered_variants/$*.indel.dbsnp.snpeff.dbNSFP.filtered.vcf \
		--rmsk-file ~/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file ~/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file ~/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file ~/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible ~/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro ~/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--rejected-variants-file ~/hdall/results/curated-recected-variants.tsv \
		--remission-variants-file remission-variants.tsv.gz \
		--evs-file ~/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		--min-num-rem-to-exclude 3 \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@	
	
#----------------
# PICARD
#----------------

.PHONY: picard
picard: 

picard/%.picard.insertsize.out: ~/p2ry8-crlf2/data/bam/%.duplicate_marked.realigned.recalibrated.bam ~/tools/picard-tools-1.114/CollectInsertSizeMetrics.jar
	java -jar ~/tools/picard-tools-1.114/CollectInsertSizeMetrics.jar \
		INPUT=$< \
		HISTOGRAM_FILE=picard/$*.picard.insertsize.pdf \
		OUTPUT=$@.part \
		STOP_AFTER=10000000
	mv $@.part $@
	
#---------------
# PINDEL
#---------------

pindel/pindel.cfg: $(foreach P, $(PATIENTS_MATCHED), picard/$PC.picard.insertsize.out picard/$PD.picard.insertsize.out picard/$PR.picard.insertsize.out) \
				   $(foreach P, $(PATIENTS_DIA_ONLY), picard/$PC.picard.insertsize.out picard/$PD.picard.insertsize.out) \
		           $(foreach P, $(PATIENTS_REL2), picard/$PR2.picard.insertsize.out) \
		           $(foreach P, $(PATIENTS_REL3), picard/$PR3.picard.insertsize.out)
	grep -A 1 MEDIAN_INSERT_SIZE $^ | perl -ne 'if (/\/([^\.]+)\.picard.insertsize.out-(\d+)/) { print "/home/STANNANET/christian.frech/p2ry8-crlf2/data/bam/$$1.duplicate_marked.realigned.recalibrated.bam\t$$2\t$$1\n"; }' > $@.part 
	mv $@.part $@ 

.PHONY: pindel
pindel: pindel/allsamples.pindel.tsv

pindel/allsamples_D pindel/allsamples_SI: pindel/pindel.cfg ~/tools/pindel-0.2.4w/pindel ~/generic/data/broad/hs37d5.fa
	~/tools/pindel-0.2.4w/pindel \
		--fasta ~/generic/data/broad/hs37d5.fa \
		--config-file pindel/pindel.cfg \
		--output-prefix pindel/allsamples \
		--report_long_insertions true \
		--NormalSamples true \
		--minimum_support_for_event 10 \
		--number_of_threads 20

pindel/allsamples.combined.filtered.vcf: pindel/allsamples_D.vcf pindel/allsamples_SI.vcf pindel/allsamples_LI.vcf pindel/allsamples_TD.vcf ~/p2ry8-crlf2/scripts/filter-pindel.pl
	~/tools/vcftools_0.1.10/bin/vcf-concat \
		<(perl ~/p2ry8-crlf2/scripts/filter-pindel.pl --vcf-in pindel/allsamples_D.vcf) \
		<(perl ~/p2ry8-crlf2/scripts/filter-pindel.pl --vcf-in pindel/allsamples_SI.vcf) \
		<(perl ~/p2ry8-crlf2/scripts/filter-pindel.pl --vcf-in pindel/allsamples_LI.vcf) \
		<(perl ~/p2ry8-crlf2/scripts/filter-pindel.pl --vcf-in pindel/allsamples_TD.vcf) \
		| ~/tools/vcftools_0.1.10/bin/vcf-sort > $@.part
	mv $@.part $@
	
pindel/%.vcf: pindel/%
	~/tools/pindel-0.2.4w/pindel2vcf -p $< -r ~/generic/data/broad/hs37d5.fa -R hs37d5.fa -d 2011-07-01 --both_strands_supported -v $@.part
	mv $@.part $@

pindel/allsamples.combined.filtered.dbsnp.vcf: pindel/allsamples.combined.filtered.vcf ~/generic/data/ncbi/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/generic/data/ncbi/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat ~/p2ry8-crlf2/results/$< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@

pindel/allsamples.combined.filtered.dbsnp.snpeff.vcf: pindel/allsamples.combined.filtered.dbsnp.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -Xmx2g -jar snpEff.jar -v -lof hg19 -stats $(PWD)/snpeff/$*.snpeff.summary.html $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

pindel/allsamples.combined.filtered.dbsnp.snpeff.dbNSFP.vcf: pindel/allsamples.combined.filtered.dbsnp.snpeff.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar dbnsfp -v ~/generic/data/dbNSFP-2.1/dbNSFP2.1.txt $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

pindel/allsamples.pindel.tsv: pindel/allsamples.combined.filtered.dbsnp.snpeff.dbNSFP.vcf ~/p2ry8-crlf2/scripts/pindel-vcf-to-tsv.pl
	perl ~/p2ry8-crlf2/scripts/pindel-vcf-to-tsv.pl \
		--vcf-in $< \
		--evs-file ~/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		> $@.part
	mv $@.part $@ 