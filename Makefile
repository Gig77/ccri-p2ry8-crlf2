export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PATIENTS_MATCHED = 108 839 92 B36 BB16 GI13 HV57 HV80 LU3 N7 S23 SN18 737
PATIENTS_DIA_ONLY = 242 360 365 379 400 506 769 833 948
PATIENTS_REL2 = 108 737 

#all: filtered-variants.cosmic.tsv filtered-variants.cosmic.normaf.tsv snpeff/mutect_737_rem_rel1.dbsnp.snpeff.dbNSFP.vcf snpeff/mutect_737_rem_rel2.dbsnp.snpeff.dbNSFP.vcf snpeff/mutect_737_rem_dia.dbsnp.snpeff.dbNSFP.vcf snpeff/indels_737_rem_rel1.indel.dbsnp.snpeff.dbNSFP.vcf snpeff/indels_737_rem_rel2.indel.dbsnp.snpeff.dbNSFP.vcf snpeff/indels_737_rem_dia.indel.dbsnp.snpeff.dbNSFP.vcf
all: filtered-variants.cosmic.normaf.tsv filtered-variants.cosmic.merged.tsv 

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
# SNPEFF
#-----------	
snpeff/%.dbsnp.vcf: ~/p2ry8-crlf2/data/mutect_somatic_mutations_hg19/%_calls.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat $< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@

snpeff/%.indel.dbsnp.vcf: ~/p2ry8-crlf2/data/mutect_somatic_indels_hg19/%_indel.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf \
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
# final lists
#-------------
filtered-variants.tsv: $(foreach P, $(PATIENTS_MATCHED), filtered_variants/$P_rem_dia.snp.filtered.tsv filtered_variants/$P_rem_dia.indel.filtered.tsv filtered_variants/$P_rem_rel.snp.filtered.tsv filtered_variants/$P_rem_rel.indel.filtered.tsv) \
					   $(foreach P, $(PATIENTS_DIA_ONLY), filtered_variants/$P_rem_dia.snp.filtered.tsv filtered_variants/$P_rem_dia.indel.filtered.tsv) \
					   $(foreach P, $(PATIENTS_REL2), filtered_variants/$P_rem_rel2.snp.filtered.tsv filtered_variants/$P_rem_rel2.indel.filtered.tsv)
#					   filtered_variants/715_rem_rel3.snp.filtered.tsv \
#					   filtered_variants/715_rem_rel3.indel.filtered.tsv \
#					   filtered_variants/737_rem_dia.indel.filtered.tsv \
#					   filtered_variants/737_rem_rel.indel.filtered.tsv \
#					   filtered_variants/737_rem_rel2.indel.filtered.tsv \
#					   filtered_variants/108_rem_dia.indel.filtered.tsv \
#					   filtered_variants/108_rem_rel.indel.filtered.tsv \
#					   filtered_variants/108_rem_rel2.indel.filtered.tsv \
#					   ~/hdall/scripts/filter-variants.pl 
	perl  ~/hdall/scripts/filter-variants.pl --header 2>&1 1>$@.part | $(LOG)
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

filtered_variants/%.snp.filtered.tsv: snpeff/mutect_%.dbsnp.snpeff.vcf ~/hdall/results/curated-recected-variants.tsv ~/hdall/results/remission-variants.tsv.gz.tbi
	mkdir -p filtered_variants
	perl ~/hdall/scripts/filter-variants.pl \
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
		--remission-variants-file ~/hdall/results/remission-variants.tsv.gz \
		--evs-file ~/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@

filtered_variants/%.indel.filtered.tsv: snpeff/indels_%.indel.dbsnp.snpeff.vcf ~/hdall/results/curated-recected-variants.tsv ~/hdall/results/remission-variants.tsv.gz.tbi
	mkdir -p filtered_variants
	perl ~/hdall/scripts/filter-variants.pl \
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
		--remission-variants-file ~/hdall/results/remission-variants.tsv.gz \
		--evs-file ~/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@	
