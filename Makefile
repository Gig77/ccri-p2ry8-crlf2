export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
.SECONDEXPANSION:

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PATIENTS_TEST = 92
PATIENTS_MATCHED = 108 92 460 545 564 715 737 839 B36 BB16 DL2 GI8 GI13 HV57 HV80 LU3 MA5 N7 S23 SN18 DS10898 VS14645 SE15285 BJ17183 KE17247 AL9890 GL11356
PATIENTS_DIA_ONLY = 242 360 365 379 400 506 769 802 833 887 841 903 948 957 961 1060 1066 1089 HW11537 KT14158 TL14516
PATIENTS_REL2 = 108 737 AL9890
PATIENTS_REL3 = 715 737 S23
PATIENTS_XENO = m1963-545-rel m1964-545-rel m1957-715-rel m1035-108-dia m1060-108-rel m252-379-dia m1041-737-dia m1069-737-rel m247-833-dia m1059-92-dia m1037-839-dia m248-841-dia m1977-G-dia m1967-Y-rel m1966-Y-dia

#all: filtered-variants.cosmic.tsv snpeff/mutect_737_rem_rel1.dbsnp.snpeff.dbNSFP.vcf snpeff/mutect_737_rem_rel2.dbsnp.snpeff.dbNSFP.vcf snpeff/mutect_737_rem_dia.dbsnp.snpeff.dbNSFP.vcf snpeff/indels_737_rem_rel1.indel.dbsnp.snpeff.dbNSFP.vcf snpeff/indels_737_rem_rel2.indel.dbsnp.snpeff.dbNSFP.vcf snpeff/indels_737_rem_dia.indel.dbsnp.snpeff.dbNSFP.vcf
all: snp fastqc gene-patient-matrix.tsv filtered-variants.paper.tsv filtered-variants.xenografts.cosmic.tsv filtered-variants.cosmic.merged.tsv coverage-region coverage-genome/allpatients.coverage-genome.pdf coverage-chr21/allpatients.coverage-chr21.pdf picard loh snp-comparison/allsamples.snp-comparison.tsv figures impacted-genes.tsv

#-----------
# DOWNLOAD
#-----------
#wget -r -np -e robots=off --reject *.wig.txt,*_stats.out http://www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_zahyee0Zooxiphoogo5ooMee3mai8uy8/mutect_b37/ 
#wget -r -np -e robots=off --reject *.wig.txt,*_stats.out http://www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_zahyee0Zooxiphoogo5ooMee3mai8uy8/indels_b37/
#MuTect,Indelocator results only:
#wget -r -np -e robots=off --accept D_annotated.vcf,R_annotated.vcf,R1_annotated.vcf,R2_annotated.vcf,Diagnosis_annotated.vcf,Relapse_annotated.vcf --limit-rate=1500k http://www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_zahyee0Zooxiphoogo5ooMee3mai8uy8/
#reprocessed HD samples
#wget -r -np -e robots=off --accept *715*,*460*,*564*,*545* --limit-rate=1500k http://www.biomedical-sequencing.at/projects/BSA_0003_HDBALL_1ca2f4fb9c73e5eb3986819f818f8acb/
#wget -r -np -e robots=off --accept *AL9890*,*GL11356*,*737_Relapse3* --limit-rate=1500k http://www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_zahyee0Zooxiphoogo5ooMee3mai8uy8/
#wget -r -np -e robots=off --accept *M1957*,*m1963*,*m1964*,*m1967*,*m1977* --limit-rate=1500k http://www.biomedical-sequencing.at/projects/BSA_0024_HD_BCP_ALL_xiewahYoojiev7ahyoh4ohjiseiKepoh/
#wget -r -np -e robots=off --accept *AL9890*,*S23*,*M248*,*M247*,*M252*,*M1059*,*M1060*,*M1069*,*M1035*,*M1037*,*M1041* --reject *.wig --limit-rate=1500k http://www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/  // PC mouse xenografts, AL9890 2nd relapse, S23 3rd relapse, downloaded Aug 26th, 2015
#wget -r -np -e robots=off --accept *M1966* --reject *.wig --limit-rate=1500k http://www.biomedical-sequencing.at/projects/BSA_0003_HDBALL_eb8f0fe088814dcd95e115a3f4ac2e74/  // mouse xeonograft from HD case Y (diagnosis), downloaded Aug 26th, 2015

/mnt/projects/p2ry8-crlf2/data/mutect_somatic_mutations_hg19/%_calls.vcf: /mnt/projects/p2ry8-crlf2/data/mutect_somatic_mutations/%_calls.vcf
	cat $< | perl -ne 's/^([\dXY])/chr$$1/; s/^MT/chrM/; print $$_;' > $@.part
	mv $@.part $@

/mnt/projects/p2ry8-crlf2/data/mutect_somatic_indels_hg19/%_indel.vcf: /mnt/projects/p2ry8-crlf2/data/mutect_somatic_indels/%_indel.vcf
	cat $< | perl -ne 's/^([\dXYM])/chr$$1/; print $$_;' > $@.part
	mv $@.part $@
	
# concatenate mutation and indel VCFs for samples sequenced in a later batch (because CeMM did not do it as for the previous batch)
#-----

/mnt/projects/p2ry8-crlf2/data/mutect/m1966-Y-dia_rem_xeno.somatic.vcf: /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0003_HDBALL_eb8f0fe088814dcd95e115a3f4ac2e74/variant_calling_somatic_Y3767_Remission__M1966_Y_Diagnosis_mutations.vcf /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0003_HDBALL_eb8f0fe088814dcd95e115a3f4ac2e74/variant_calling_somatic_Y3767_Remission__M1966_Y_Diagnosis_indels.vcf
	(head -1 $(word 1, $^); grep ^## $(word 1, $^) | tail -n +2; grep ^## $(word 2, $^) | tail -n +2; grep ^#CHROM $(word 1, $^); grep -v ^# $(word 1, $^); grep -v ^# $(word 2, $^)) > $@
/mnt/projects/p2ry8-crlf2/data/mutect/m1035-108-dia_rem_xeno.somatic.vcf: /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_108C__M1035_108_Diagnosis_mutations.vcf /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_108C__M1035_108_Diagnosis_indels.vcf
	(head -1 $(word 1, $^); grep ^## $(word 1, $^) | tail -n +2; grep ^## $(word 2, $^) | tail -n +2; grep ^#CHROM $(word 1, $^); grep -v ^# $(word 1, $^); grep -v ^# $(word 2, $^)) > $@
/mnt/projects/p2ry8-crlf2/data/mutect/m1060-108-rel_rem_xeno.somatic.vcf: /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_108C__M1060_108_Relapse_mutations.vcf /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_108C__M1060_108_Relapse_indels.vcf
	(head -1 $(word 1, $^); grep ^## $(word 1, $^) | tail -n +2; grep ^## $(word 2, $^) | tail -n +2; grep ^#CHROM $(word 1, $^); grep -v ^# $(word 1, $^); grep -v ^# $(word 2, $^)) > $@
/mnt/projects/p2ry8-crlf2/data/mutect/m252-379-dia_rem_xeno.somatic.vcf: /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_379C__M252_379_Diagnosis_mutations.vcf /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_379C__M252_379_Diagnosis_indels.vcf
	(head -1 $(word 1, $^); grep ^## $(word 1, $^) | tail -n +2; grep ^## $(word 2, $^) | tail -n +2; grep ^#CHROM $(word 1, $^); grep -v ^# $(word 1, $^); grep -v ^# $(word 2, $^)) > $@
/mnt/projects/p2ry8-crlf2/data/mutect/m1041-737-dia_rem_xeno.somatic.vcf: /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_737C__M1041_737_Diagnosis_mutations.vcf /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_737C__M1041_737_Diagnosis_indels.vcf
	(head -1 $(word 1, $^); grep ^## $(word 1, $^) | tail -n +2; grep ^## $(word 2, $^) | tail -n +2; grep ^#CHROM $(word 1, $^); grep -v ^# $(word 1, $^); grep -v ^# $(word 2, $^)) > $@
/mnt/projects/p2ry8-crlf2/data/mutect/m1069-737-rel_rem_xeno.somatic.vcf: /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_737C__M1069_737_Relapse_mutations.vcf /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_737C__M1069_737_Relapse_indels.vcf
	(head -1 $(word 1, $^); grep ^## $(word 1, $^) | tail -n +2; grep ^## $(word 2, $^) | tail -n +2; grep ^#CHROM $(word 1, $^); grep -v ^# $(word 1, $^); grep -v ^# $(word 2, $^)) > $@
/mnt/projects/p2ry8-crlf2/data/mutect/m247-833-dia_rem_xeno.somatic.vcf: /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_833C__M247_833_Diagnosis_mutations.vcf /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_833C__M247_833_Diagnosis_indels.vcf
	(head -1 $(word 1, $^); grep ^## $(word 1, $^) | tail -n +2; grep ^## $(word 2, $^) | tail -n +2; grep ^#CHROM $(word 1, $^); grep -v ^# $(word 1, $^); grep -v ^# $(word 2, $^)) > $@
/mnt/projects/p2ry8-crlf2/data/mutect/m1037-839-dia_rem_xeno.somatic.vcf: /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_839C__M1037_839_Diagnosis_mutations.vcf /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_839C__M1037_839_Diagnosis_indels.vcf
	(head -1 $(word 1, $^); grep ^## $(word 1, $^) | tail -n +2; grep ^## $(word 2, $^) | tail -n +2; grep ^#CHROM $(word 1, $^); grep -v ^# $(word 1, $^); grep -v ^# $(word 2, $^)) > $@
/mnt/projects/p2ry8-crlf2/data/mutect/m248-841-dia_rem_xeno.somatic.vcf: /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_842_Remission__M248_841_Diagnosis_mutations.vcf /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_842_Remission__M248_841_Diagnosis_indels.vcf
	(head -1 $(word 1, $^); grep ^## $(word 1, $^) | tail -n +2; grep ^## $(word 2, $^) | tail -n +2; grep ^#CHROM $(word 1, $^); grep -v ^# $(word 1, $^); grep -v ^# $(word 2, $^)) > $@
/mnt/projects/p2ry8-crlf2/data/mutect/m1059-92-dia_rem_xeno.somatic.vcf: /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_92C__M1059_92_Diagnosis_mutations.vcf /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_92C__M1059_92_Diagnosis_indels.vcf
	(head -1 $(word 1, $^); grep ^## $(word 1, $^) | tail -n +2; grep ^## $(word 2, $^) | tail -n +2; grep ^#CHROM $(word 1, $^); grep -v ^# $(word 1, $^); grep -v ^# $(word 2, $^)) > $@
/mnt/projects/p2ry8-crlf2/data/mutect/AL9890_rem_rel2.somatic.vcf: /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_AL9890_Remission__AL9890_RR_mutations.vcf /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_AL9890_Remission__AL9890_RR_indels.vcf
	(head -1 $(word 1, $^); grep ^## $(word 1, $^) | tail -n +2; grep ^## $(word 2, $^) | tail -n +2; grep ^#CHROM $(word 1, $^); grep -v ^# $(word 1, $^); grep -v ^# $(word 2, $^)) > $@
/mnt/projects/p2ry8-crlf2/data/mutect/S23_rem_rel3.somatic.vcf: /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_S23C__S23_RR_mutations.vcf /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_b4e692fc2393447aadea66b2bfeb2863/variant_calling_somatic_S23C__S23_RR_indels.vcf
	(head -1 $(word 1, $^); grep ^## $(word 1, $^) | tail -n +2; grep ^## $(word 2, $^) | tail -n +2; grep ^#CHROM $(word 1, $^); grep -v ^# $(word 1, $^); grep -v ^# $(word 2, $^)) > $@
	
#-----------
# REMISSION VARIANTS
#-----------
remission-variants.tsv: $(foreach P, $(PATIENTS_MATCHED), remission-variants/$P.remission-variants.tsv) \
					    $(foreach P, $(PATIENTS_DIA_ONLY), remission-variants/$P.remission-variants.tsv) \
					    /mnt/projects/hdall/results/remission-variants.tsv
	cat $^ | perl -ne 's/chrM/MT/; s/\tchr/\t/; print $$_;' | ~/tools/lh3-sort/sort -k 2,2N -k 3,3n > $@.part
	mv $@.part $@

remission-variants.tsv.gz: remission-variants.tsv
	bgzip -c $^ > $@.part
	mv $@.part $@

remission-variants.tsv.gz.tbi: remission-variants.tsv.gz
	/data_synology/software/tabix-0.2.6/tabix $^ -s 2 -b 3 -e 3

remission-variants/%.remission-variants.tsv: /mnt/projects/p2ry8-crlf2/data/mutect/%_rem_dia.somatic.vcf /mnt/projects/p2ry8-crlf2/scripts/create-normal-vcf.pl
	mkdir -p remission-variants 
	cat $< | perl /mnt/projects/p2ry8-crlf2/scripts/create-normal-vcf.pl --sample-name $* \
		2>&1 1> $@.part | $(LOG)
	mv $@.part $@

#-----------	
# SNPEFF
#-----------	
snpeff/%.dbsnp.vcf: /mnt/projects/p2ry8-crlf2/data/mutect/%.vcf /mnt/projects/generic/data/ncbi/common_no_known_medical_impact_20140826.vcf
	PWD=$(pwd)
	(cd /data_synology/software/snpeff-3.6; java -jar SnpSift.jar annotate \
		-v /mnt/projects/generic/data/ncbi/common_no_known_medical_impact_20140826.vcf \
		<(cat $< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@

snpeff/%.dbsnp.snpeff.vcf: snpeff/%.dbsnp.vcf
	PWD=$(pwd)
	(cd /data_synology/software/snpeff-3.6; java -Xmx4g -jar snpEff.jar -v -lof GRCh37.75 -stats $(PWD)/snpeff/$*.snpeff.summary.html $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

snpeff/%.dbsnp.snpeff.dbNSFP.vcf: snpeff/%.dbsnp.snpeff.vcf /mnt/projects/generic/data/dbNSFP-2.6/dbNSFP2.6_variant.tsv.gz
	PWD=$(pwd)
	(cd /data_synology/software/snpeff-3.6; java -jar SnpSift.jar dbnsfp \
		-v /mnt/projects/generic/data/dbNSFP-2.6/dbNSFP2.6_variant.tsv.gz \
		-collapse \
		-f SIFT_pred,SIFT_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_score,SiPhy_29way_logOdds,LRT_pred,LRT_score,MutationTaster_pred,MutationTaster_score,MutationAssessor_pred,MutationAssessor_score,FATHMM_pred,FATHMM_score,RadialSVM_pred,RadialSVM_score,GERP++_RS,1000Gp1_AF,1000Gp1_AFR_AF,1000Gp1_EUR_AF,1000Gp1_AMR_AF,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF,Uniprot_acc,Interpro_domain, \
		$(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

#-------------
# FASTQC
#-------------

.PHONY: fastqc
fastqc: $(foreach P, $(PATIENTS_MATCHED), fastqc/variant_calling_process_sample_$PC_realigned_fastqc.html fastqc/variant_calling_process_sample_$PD_realigned_fastqc.html fastqc/variant_calling_process_sample_$PR_realigned_fastqc.html) \
        $(foreach P, $(PATIENTS_DIA_ONLY), fastqc/variant_calling_process_sample_$PC_realigned_fastqc.html fastqc/variant_calling_process_sample_$PD_realigned_fastqc.html) \
		$(foreach P, $(PATIENTS_REL2), fastqc/variant_calling_process_sample_$PR2_realigned_fastqc.html) \
		$(foreach P, $(PATIENTS_REL3), fastqc/variant_calling_process_sample_$PR3_realigned_fastqc.html) \
		$(foreach P, $(PATIENTS_XENO), fastqc/variant_calling_process_sample_$P_realigned_fastqc.html)
	
fastqc/variant_calling_process_sample_%_realigned_fastqc.html: /mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_%_realigned.bam
	mkdir -p fastqc/$*.part
	/data_synology/software/FastQC-0.11.2/fastqc -o fastqc/$*.part -f bam $^
	mv fastqc/$*.part/* fastqc
	rmdir fastqc/$*.part

#-------------
# REGION COVERAGE
#-------------

GENES_PAR1=PPP2R3B,SHOX,CRLF2,IL3RA,P2RY8,ASMT,DHR3X,ZBED1,CD99
GENES_IKZF1=VWC2,ZPBP,C7orf72,IKZF1,DDC,GRB10,COBL

.PHONY: coverage-region
coverage-region: coverage-region/allpatients.PAR1-X-60001-2699520.pdf \
		  coverage-region/allpatients.IKZF1-7-49578046-51601231.pdf \
		  coverage-region/allpatients.IKZF1_highres-7-50332565-50494236.pdf \
		  coverage-region/allpatients.IKZF2-2-211852462-217849831.pdf \
		  coverage-region/allpatients.IKZF3-17-37789665-38184756.pdf \
		  coverage-region/allpatients.PAX5-9-35800105-38084658.pdf \
		  coverage-region/allpatients.EBF1-5-155175854-161576670.pdf \
		  coverage-region/allpatients.ETV6-12-10309186-13646908.pdf \
		  coverage-region/allpatients.RUNX1-21-34971197-38263248.pdf \
		  coverage-region/allpatients.VPREB1-22-21866719-24317257.pdf \
		  coverage-region/allpatients.ERG-21-38349519-41562975.pdf \
		  coverage-region/allpatients.ERG_highres-21-39739183-40033618.pdf \
		  coverage-region/allpatients.TP53-17-7440139-7721205.pdf \
		  coverage-region/allpatients.RB1-13-46907242-51051394.pdf \
		  coverage-region/allpatients.CDKN2AandB-9-20309360-23721195.pdf \
		  coverage-region/allpatients.CREBBP-16-3067833-4618486.pdf \
		  coverage-region/allpatients.MLL2-12-49295653-49584389.pdf \
		  coverage-region/allpatients.EZH2-7-148002918-149209618.pdf \
		  coverage-region/allpatients.NCOR1-17-15245025-16768909.pdf \
		  coverage-region/allpatients.TUSC3-8-12275019-19439531.pdf \
		  coverage-region/allpatients.WHSC1-4-1615018-2430195.pdf \
		  coverage-region/allpatients.NT5C2-10-104141780-105714035.pdf \
		  coverage-region/allpatients.LEF1-4-108313628-110095704.pdf \
		  coverage-region/allpatients.TCF3-19-1261454-2018719.pdf \
		  coverage-region/allpatients.BLNK-10-97442705-98708695.pdf \
		  coverage-region/allpatients.FOXO3A-6-108206749-109779719.pdf \
		  coverage-region/allpatients.FBXW7-4-151992600-154510682.pdf \
		  coverage-region/allpatients.CREG1-1-167207942-167849758.pdf \
		  coverage-region/allpatients.FLT3_PAN3-13-27955721-29384162.pdf \
		  coverage-region/allpatients.PDGFRB-5-149383560-149648693.pdf \
		  coverage-region/allpatients.STRN3-14-30942316-31995615.pdf \
		  coverage-region/allpatients.RANBP2-2-109007910-109671216.pdf \
		  coverage-region/allpatients.EPOR-19-11393474-11619196.pdf \
		  coverage-region/allpatients.SH2B3-12-111621200-112196224.pdf \
		  coverage-region/allpatients.HIST1H2BD_HIST1H1E-6-26033646-26316674.pdf \
		  coverage-region/allpatients.SPRED1-15-35957458-40803200.pdf \
		  coverage-region/allpatients.ADD3-10-111062840-113085360.pdf \
		  coverage-region/allpatients.ATP10A-15-24452756-27861838.pdf \
		  coverage-region/allpatients.NUP214-9-133609878-134537138.pdf \
		  coverage-region/allpatients.BTG1-12-91300000-93550000.pdf \
		  coverage-region/allpatients.GABRB3-15-25870574-27936343.pdf \
		  coverage-region/allpatients.SETD2-3-46897168-47324570.pdf
		  
coverage-region/%.exon-coverage.tsv: /mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_%_realigned.bam /mnt/projects/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed /mnt/projects/hdall/scripts/cnv/get-exon-coverage.pl
	/data_synology/software/samtools-0.1.19/samtools depth \
		-Q 1 \
		-b /mnt/projects/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed \
		$< \
	| perl ~/git/hdall/cnv/get-exon-coverage.pl \
		--exon-bed /mnt/projects/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@

coverage-region/patient%.matched.pdf: coverage-region/$$(word 1, $$(subst ., , %))D.exon-coverage.tsv coverage-region/$$(word 1, $$(subst ., , %))R.exon-coverage.tsv coverage-region/$$(word 1, $$(subst ., , %))C.exon-coverage.tsv /mnt/projects/p2ry8-crlf2/scripts/cov-plot-region.R
	Rscript /mnt/projects/p2ry8-crlf2/scripts/cov-plot-region.R \
		--patient $(word 1, $(subst ., , $*)) \
		--sample1 $(word 1,$^) \
		--sample2 $(word 2,$^) \
		--remission $(word 3,$^) \
		--output $@.part \
		--region-name $(word 1, $(subst -, , $(word 2, $(subst ., , $*)))) \
		--display-chrom $(word 2, $(subst -, , $*)) \
		--display-start $(word 3, $(subst -, , $*)) \
		--display-end $(word 4, $(subst -, , $*)) \
		$(if $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),--display-genes $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),)
	mv $@.part $@

coverage-region/patient%.diaonly.pdf: coverage-region/$$(word 1, $$(subst ., , %))D.exon-coverage.tsv coverage-region/$$(word 1, $$(subst ., , %))C.exon-coverage.tsv /mnt/projects/p2ry8-crlf2/scripts/cov-plot-region.R
	Rscript /mnt/projects/p2ry8-crlf2/scripts/cov-plot-region.R \
		--patient $(word 1, $(subst ., , $*)) \
		--sample1 $(word 1,$^) \
		--remission $(word 2,$^) \
		--output $@.part \
		--region-name $(word 1, $(subst -, , $(word 2, $(subst ., , $*)))) \
		--display-chrom $(word 2, $(subst -, , $*)) \
		--display-start $(word 3, $(subst -, , $*)) \
		--display-end $(word 4, $(subst -, , $*)) \
		$(if $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),--display-genes $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),)
	mv $@.part $@

coverage-region/patient%.rel2.pdf: coverage-region/$$(word 1, $$(subst ., , %))R2.exon-coverage.tsv coverage-region/$$(word 1, $$(subst ., , %))C.exon-coverage.tsv /mnt/projects/p2ry8-crlf2/scripts/cov-plot-region.R
	Rscript /mnt/projects/p2ry8-crlf2/scripts/cov-plot-region.R \
		--patient $(word 1, $(subst ., , $*)) \
		--sample1 $(word 1,$^) \
		--remission $(word 2,$^) \
		--output $@.part \
		--region-name $(word 1, $(subst -, , $(word 2, $(subst ., , $*)))) \
		--display-chrom $(word 2, $(subst -, , $*)) \
		--display-start $(word 3, $(subst -, , $*)) \
		--display-end $(word 4, $(subst -, , $*)) \
		$(if $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),--display-genes $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),)
	mv $@.part $@

coverage-region/patient%.rel3.pdf: coverage-region/$$(word 1, $$(subst ., , %))R3.exon-coverage.tsv coverage-region/$$(word 1, $$(subst ., , %))C.exon-coverage.tsv /mnt/projects/p2ry8-crlf2/scripts/cov-plot-region.R
	Rscript /mnt/projects/p2ry8-crlf2/scripts/cov-plot-region.R \
		--patient $(word 1, $(subst ., , $*)) \
		--sample1 $(word 1,$^) \
		--remission $(word 2,$^) \
		--output $@.part \
		--region-name $(word 1, $(subst -, , $(word 2, $(subst ., , $*)))) \
		--display-chrom $(word 2, $(subst -, , $*)) \
		--display-start $(word 3, $(subst -, , $*)) \
		--display-end $(word 4, $(subst -, , $*)) \
		$(if $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),--display-genes $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),)
	mv $@.part $@

coverage-region/allpatients.%.pdf: $(foreach P, $(PATIENTS_MATCHED), coverage-region/patient$P.%.matched.pdf) $(foreach P, $(PATIENTS_DIA_ONLY), coverage-region/patient$P.%.diaonly.pdf) $(foreach P, $(PATIENTS_REL2), coverage-region/patient$P.%.rel2.pdf) $(foreach P, $(PATIENTS_REL3), coverage-region/patient$P.%.rel3.pdf)
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@
	rm $^

#-------------
# WHOLE-GENOME COVERAGE
#-------------

coverage-genome/allpatients.coverage-genome.pdf: $(foreach P, $(PATIENTS_MATCHED), coverage-genome/$PD.coverage-genome.pdf coverage-genome/$PR.coverage-genome.pdf) \
								    		     $(foreach P, $(PATIENTS_DIA_ONLY), coverage-genome/$PD.coverage-genome.pdf) \
								    		     $(foreach P, $(PATIENTS_REL2), coverage-genome/$PR2.coverage-genome.pdf) \
								    		     $(foreach P, $(PATIENTS_REL3), coverage-genome/$PR3.coverage-genome.pdf) \
								    		     $(foreach P, $(PATIENTS_XENO), coverage-genome/$P.xeno.coverage-genome.pdf)
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@
	rm $^

coverage-genome/%D.coverage-genome.pdf: coverage-genome/%D.coverage-genome.tsv coverage-genome/%C.coverage-genome.tsv /mnt/projects/hdall/scripts/cnv/plot-segment-coverage.R
	mkdir -p coverage-genome/circos
	Rscript /mnt/projects/hdall/scripts/cnv/plot-segment-coverage.R \
		--patient $*D \
		--tumor $(word 1,$^) \
		--normal $(word 2,$^) \
		--circos coverage-genome/circos/$*D.cnv.circos.tsv \
		--output $@.part
	mv $@.part $@

coverage-genome/%R.coverage-genome.pdf: coverage-genome/%R.coverage-genome.tsv coverage-genome/%C.coverage-genome.tsv /mnt/projects/hdall/scripts/cnv/plot-segment-coverage.R
	mkdir -p coverage-genome/circos
	Rscript /mnt/projects/hdall/scripts/cnv/plot-segment-coverage.R \
		--patient $*R \
		--tumor $(word 1,$^) \
		--normal $(word 2,$^) \
		--circos coverage-genome/circos/$*R.cnv.circos.tsv \
		--output $@.part
	mv $@.part $@

coverage-genome/%R2.coverage-genome.pdf: coverage-genome/%R2.coverage-genome.tsv coverage-genome/%C.coverage-genome.tsv /mnt/projects/hdall/scripts/cnv/plot-segment-coverage.R
	mkdir -p coverage-genome/circos
	Rscript /mnt/projects/hdall/scripts/cnv/plot-segment-coverage.R \
		--patient $*R2 \
		--tumor $(word 1,$^) \
		--normal $(word 2,$^) \
		--circos coverage-genome/circos/$*R2.cnv.circos.tsv \
		--output $@.part
	mv $@.part $@

coverage-genome/%R3.coverage-genome.pdf: coverage-genome/%R3.coverage-genome.tsv coverage-genome/%C.coverage-genome.tsv /mnt/projects/hdall/scripts/cnv/plot-segment-coverage.R
	mkdir -p coverage-genome/circos
	Rscript /mnt/projects/hdall/scripts/cnv/plot-segment-coverage.R \
		--patient $*R3 \
		--tumor $(word 1,$^) \
		--normal $(word 2,$^) \
		--circos coverage-genome/circos/$*R3.cnv.circos.tsv \
		--output $@.part
	mv $@.part $@

coverage-genome/%.xeno.coverage-genome.pdf: coverage-genome/%.coverage-genome.tsv coverage-genome/948C.coverage-genome.tsv
	mkdir -p coverage-genome/circos
	Rscript /mnt/projects/hdall/scripts/cnv/plot-segment-coverage.R \
		--patient $* \
		--tumor $(word 1,$^) \
		--normal $(word 2,$^) \
		--circos coverage-genome/circos/$*.cnv.circos.tsv \
		--output $@.part
	mv $@.part $@

coverage-genome/%.coverage-genome.tsv: /mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_%_realigned.bam
	mkdir -p coverage-genome
	/data_synology/software/samtools-0.1.19/samtools depth -Q 1 $< \
		| perl /mnt/projects/hdall/scripts/cnv/get-segment-coverage.pl --sample $* --bin-size 250000 --chr-sizes /mnt/projects/generic/data/hg19/ucsc.hg19.chrom.sizes --add-chr \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@

#-------------
# COVERAGE CHROMOSOME 21
#-------------

coverage-chr21/allpatients.coverage-chr21.pdf: $(foreach P, $(PATIENTS_MATCHED) $(PATIENTS_DIA_ONLY), coverage-chr21/$PC.coverage-chr21.pdf)
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@
	rm $^

coverage-chr21/%C.coverage-chr21.pdf: coverage-genome/%C.coverage-genome.tsv coverage-genome/B36C.coverage-genome.tsv /mnt/projects/hdall/scripts/cnv/plot-segment-coverage-chr21.R
	mkdir -p coverage-chr21 
	Rscript /mnt/projects/hdall/scripts/cnv/plot-segment-coverage-chr21.R \
		--patient $*C \
		--case $(word 1,$^) \
		--control coverage-genome/B36C.coverage-genome.tsv \
		--output $@.part
	mv $@.part $@

#-------------
# final lists
#-------------
filtered-variants.tsv: $(foreach P, $(PATIENTS_MATCHED), filtered-variants/$P_rem_dia.filtered.tsv filtered-variants/$P_rem_rel.filtered.tsv) \
					   $(foreach P, $(PATIENTS_DIA_ONLY), filtered-variants/$P_rem_dia.filtered.tsv) \
					   $(foreach P, $(PATIENTS_REL2), filtered-variants/$P_rem_rel2.filtered.tsv) \
					   $(foreach P, $(PATIENTS_REL3), filtered-variants/$P_rem_rel3.filtered.tsv) \
					   /mnt/projects/p2ry8-crlf2/scripts/filter-variants.pl 
	perl /mnt/projects/p2ry8-crlf2/scripts/filter-variants.pl --header 2>&1 1>$@.part | $(LOG)
	cat $(filter-out $(lastword $^), $^) >> $@.part
	mv $@.part $@

filtered-variants.paper.tsv: filtered-variants.tsv /mnt/projects/p2ry8-crlf2/scripts/filter-variants-paper.R 
	Rscript /mnt/projects/p2ry8-crlf2/scripts/filter-variants-paper.R 
	mv $@.part $@
	
filtered-variants.xenografts.tsv: $(foreach P, $(PATIENTS_XENO), filtered-variants/$P_rem_xeno.filtered.tsv) /mnt/projects/p2ry8-crlf2/scripts/filter-variants.pl 
	perl /mnt/projects/p2ry8-crlf2/scripts/filter-variants.pl --header 2>&1 1>$@.part | $(LOG)
	cat $(filter-out $(lastword $^), $^) >> $@.part
	mv $@.part $@

%.cosmic.tsv: %.tsv /mnt/projects/generic/data/cosmic/v73/CosmicMutantExport.tsv /mnt/projects/p2ry8-crlf2/scripts/annotate-cosmic.pl
	cat $(word 1,$^) | perl /mnt/projects/p2ry8-crlf2/scripts/annotate-cosmic.pl \
		--cosmic-mutation-file $(word 2,$^) \
		--only-confirmed \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@ 
	
filtered-variants.cosmic.merged.tsv: filtered-variants.cosmic.tsv filtered-variants.xenografts.tsv
	Rscript /mnt/projects/p2ry8-crlf2/scripts/merge-filtered-variants.R 2>&1 | $(LOG)

filtered-variants/%.filtered.tsv: snpeff/%.somatic.dbsnp.snpeff.dbNSFP.vcf /mnt/projects/hdall/results/curated-recected-variants.tsv remission-variants.tsv.gz.tbi /mnt/projects/p2ry8-crlf2/scripts/filter-variants.pl
	mkdir -p filtered-variants
	perl /mnt/projects/p2ry8-crlf2/scripts/filter-variants.pl \
		--sample $* \
		--vcf-in $< \
		--vcf-out filtered-variants/$*.filtered.vcf \
		--rmsk-file /mnt/projects/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file /mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file /mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file /mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro /mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--rejected-variants-file /mnt/projects/hdall/results/curated-recected-variants.tsv \
		--remission-variants-file remission-variants.tsv.gz \
		--evs-file /mnt/projects/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		--min-num-rem-to-exclude 3 \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@

impacted-genes.tsv: filtered-variants.tsv /mnt/projects/p2ry8-crlf2/scripts/impacted-genes.pl
	cat $< | perl /mnt/projects/p2ry8-crlf2/scripts/impacted-genes.pl > $@.part
	mv $@.part $@

gene-patient-matrix.tsv: filtered-variants.tsv /mnt/projects/p2ry8-crlf2/scripts/get-gene-patient-matrix.R
	Rscript /mnt/projects/p2ry8-crlf2/scripts/get-gene-patient-matrix.R
	mv $@.part $@

#----------------
# FIGURES
#----------------

.PHONY: figures
figures: figures/mutations-per-patient.pdf figures/cnv-plot-relapsing.pdf figures/cnv-plot-nonrelapsing.pdf figures/coverage.png figures/signaling-barplot.non-relapsing.pdf figures/signaling-barplot.matched.pdf figures/af-plot.pdf

figures/mutations-per-patient.pdf: filtered-variants.tsv /mnt/projects/p2ry8-crlf2/scripts/figures/mutations-per-patient.R
	Rscript /mnt/projects/p2ry8-crlf2/scripts/figures/mutations-per-patient.R

figures/cnv-plot-relapsing.pdf figures/cnv-plot-nonrelapsing.pdf: exomeCopy/allpatients.filtered-segments.exomeCopy.tsv /mnt/projects/p2ry8-crlf2/scripts/figures/cnv-plot-genome.R
	Rscript /mnt/projects/p2ry8-crlf2/scripts/figures/cnv-plot-genome.R

figures/signaling-barplot.non-relapsing.pdf figures/signaling-barplot.matched.pdf: hotspot-mutations.tsv filtered-variants.cosmic.tsv /mnt/projects/p2ry8-crlf2/scripts/figures/signaling-barplot.R
	Rscript /mnt/projects/p2ry8-crlf2/scripts/figures/signaling-barplot.R

figures/af-plot.pdf: filtered-variants.tsv /mnt/projects/p2ry8-crlf2/scripts/figures/af-plot.R
	Rscript /mnt/projects/p2ry8-crlf2/scripts/figures/af-plot.R
	
figures/coverage.png: $(foreach P, $(PATIENTS_MATCHED), exome-coverage/$PD.coverage.bedtools.txt exome-coverage/$PR.coverage.bedtools.txt exome-coverage/$PC.coverage.bedtools.txt) \
               					   $(foreach P, $(PATIENTS_DIA_ONLY), exome-coverage/$PD.coverage.bedtools.txt exome-coverage/$PC.coverage.bedtools.txt) \
								   $(foreach P, $(PATIENTS_REL2), exome-coverage/$PR2.coverage.bedtools.txt) \
								   $(foreach P, $(PATIENTS_REL3), exome-coverage/$PR3.coverage.bedtools.txt) \
               					   /mnt/projects/p2ry8-crlf2/scripts/figures/coverage-plot.R
	Rscript /mnt/projects/p2ry8-crlf2/scripts/figures/coverage-plot.R

exome-coverage/%.coverage.bedtools.txt: /mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_%_realigned.bam /mnt/projects/generic/data/illumina/nexterarapidcapture_exome_targetedregions.nochr.bed
	/data_synology/software/samtools-0.1.19/samtools view -bq 1 -F 0x400 $< | /data_synology/software/bedtools-2.17.0/bin/bedtools coverage -hist -abam - -b $(word 2, $^) | grep ^all > $@.part
	mv $@.part $@
		
#----------------
# PICARD
#----------------

.PHONY: picard
picard: $(foreach P, $(PATIENTS_MATCHED), picard/$PC.picard.insertsize.out picard/$PD.picard.insertsize.out picard/$PR.picard.insertsize.out) \
        $(foreach P, $(PATIENTS_DIA_ONLY), picard/$PC.picard.insertsize.out picard/$PD.picard.insertsize.out) \
		$(foreach P, $(PATIENTS_REL2), picard/$PR2.picard.insertsize.out) \
		$(foreach P, $(PATIENTS_REL3), picard/$PR3.picard.insertsize.out) \
		$(foreach P, $(PATIENTS_XENO), picard/$P.picard.insertsize.out) \
		picard/second_batch.insertsize.pdf \
		picard/allsamples.hs_metrics.picard.tsv
	
picard/%.picard.insertsize.out picard/%.picard.insertsize.pdf: /mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_%_realigned.bam /data_synology/software/picard-tools-1.114/CollectInsertSizeMetrics.jar
	mkdir -p picard
	java -jar /data_synology/software/picard-tools-1.114/CollectInsertSizeMetrics.jar \
		INPUT=$< \
		HISTOGRAM_FILE=picard/$*.picard.insertsize.pdf.part \
		OUTPUT=picard/$*.picard.insertsize.out.part \
		STOP_AFTER=10000000
	mv picard/$*.picard.insertsize.pdf.part picard/$*.picard.insertsize.pdf
	mv picard/$*.picard.insertsize.out.part picard/$*.picard.insertsize.out

picard/%.hs_metrics: /mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_%_realigned.bam /mnt/projects/generic/data/illumina/nexterarapidcapture_exome_targetedregions.nochr.bed
	mkdir -p picard
	/data_synology/software/samtools-0.1.19/samtools view -H $< 2>&1 1> picard/$*.targetedregions.for-picard.bed | $(LOG)
	gawk 'BEGIN { OFS="\t"} {print $$1,$$2,$$3,"+",$$4 }' /mnt/projects/generic/data/illumina/nexterarapidcapture_exome_targetedregions.nochr.bed >> picard/$*.targetedregions.for-picard.bed
	java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /data_synology/software/picard-tools-1.114/CalculateHsMetrics.jar \
		BAIT_INTERVALS=picard/$*.targetedregions.for-picard.bed \
		TARGET_INTERVALS=picard/$*.targetedregions.for-picard.bed \
		INPUT=$< \
		OUTPUT=$@.part \
		REFERENCE_SEQUENCE=/mnt/projects/generic/data/broad/hs37d5.fa \
		PER_TARGET_COVERAGE=picard/$*.hs_metrics.per_target_coverage.part \
		VALIDATION_STRINGENCY=LENIENT \
		2>&1 | $(LOG)
	mv picard/$*.hs_metrics.per_target_coverage.part picard/$*.hs_metrics.per_target_coverage
	rm picard/$*.targetedregions.for-picard.bed
	mv $@.part $@

picard/allsamples.hs_metrics.picard.tsv:  $(foreach P, $(PATIENTS_MATCHED), picard/$PC.hs_metrics picard/$PD.hs_metrics picard/$PR.hs_metrics) \
       									  $(foreach P, $(PATIENTS_DIA_ONLY), picard/$PC.hs_metrics picard/$PD.hs_metrics) \
       									  $(foreach P, $(PATIENTS_REL2), picard/$PR2.hs_metrics) \
       									  $(foreach P, $(PATIENTS_REL3), picard/$PR3.hs_metrics) \
       									  $(foreach P, $(PATIENTS_XENO), picard/$P.hs_metrics) 
	grep BAIT_SET picard/108D.hs_metrics > $@.part
	grep 3137454505 $^ | cut -f 2 -d ':' >> $@.part
	mv $@.part $@
	
picard/second_batch.insertsize.pdf: $(foreach P, m1966-Y-dia m1035-108-dia m1060-108-rel m252-379-dia m1041-737-dia m1069-737-rel m247-833-dia m1059-92-dia m1037-839-dia m248-841-dia AL9890R2 S23R3, picard/$P.picard.insertsize.pdf)
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@
	rm $^

#---------------
# PINDEL
#---------------

.PHONY: pindel
pindel: $(foreach C, 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y, pindel/chr$C_D)

pindel/pindel.cfg: $(foreach P, $(PATIENTS_MATCHED), picard/$PC.picard.insertsize.out picard/$PD.picard.insertsize.out picard/$PR.picard.insertsize.out) \
				   $(foreach P, $(PATIENTS_DIA_ONLY), picard/$PC.picard.insertsize.out picard/$PD.picard.insertsize.out) \
		           $(foreach P, $(PATIENTS_REL2), picard/$PR2.picard.insertsize.out) \
		           $(foreach P, $(PATIENTS_REL3), picard/$PR3.picard.insertsize.out)
	grep -A 1 MEDIAN_INSERT_SIZE $^ | perl -ne 'if (/\/([^\.]+)\.picard.insertsize.out-(\d+)/) { print "/mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_$$1_realigned.bam\t$$2\t$$1\n"; }' > $@.part 
	mv $@.part $@ 

#pindel/pindel.cfg: $(foreach P, 108, picard/$PC.picard.insertsize.out picard/$PD.picard.insertsize.out picard/$PR.picard.insertsize.out)
#	grep -A 1 MEDIAN_INSERT_SIZE $^ | perl -ne 'if (/\/([^\.]+)\.picard.insertsize.out-(\d+)/) { print "/mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_$$1_realigned.bam\t$$2\t$$1\n"; }' > $@.part 
#	mv $@.part $@ 

pindel/chr%_D: pindel/pindel.cfg /data_synology/software/pindel-0.2.4w/pindel /mnt/projects/generic/data/broad/hs37d5.fa
	/data_synology/software/pindel-0.2.4w/pindel \
		--fasta /mnt/projects/generic/data/broad/hs37d5.fa \
		--config-file pindel/pindel.cfg \
		--chromosome $* \
		--output-prefix pindel/chr$* \
		--report_long_insertions true \
		--NormalSamples true \
		--minimum_support_for_event 10 \
		--number_of_threads 10 \
		--max_range_index 4 \
		--report_interchromosomal_events false \
		--report_duplications false \
		--report_inversions false \
		--exclude /mnt/projects/p2ry8-crlf2/scripts/pindel-exclude-regions.bed \
		--window_size 3 \
		2>&1 | $(LOG)

# HACK; see http://stackoverflow.com/questions/19822435/multiple-targets-from-one-recipe-and-parallel-execution why this two-step process is necessary for explicit multiple target rules
pindel/allsamples_SI pindel/allsamples_LI pindel/allsamples_TD: pindel/allsamples_D

pindel/allsamples.combined.filtered.vcf: $(foreach C, 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y, pindel/chr$C_D.vcf pindel/chr$C_SI.vcf) /mnt/projects/p2ry8-crlf2/scripts/filter-pindel.pl
	/data_synology/software/vcftools_0.1.10/bin/vcf-concat \
		$(foreach C, 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y, pindel/chr$C_D.vcf) \
		$(foreach C, 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y, pindel/chr$C_SI.vcf) \
	| perl /mnt/projects/p2ry8-crlf2/scripts/filter-pindel.pl --vcf-in /dev/stdin \
	| /data_synology/software/vcftools_0.1.10/bin/vcf-sort > $@.part
	mv $@.part $@

pindel/%.vcf: pindel/%
	/data_synology/software/pindel-0.2.4w/pindel2vcf -p $< -r /mnt/projects/generic/data/broad/hs37d5.fa -R hs37d5.fa -d 2011-07-01 --both_strands_supported -v $@.part
	mv $@.part $@

pindel/allsamples.combined.filtered.dbsnp.vcf: pindel/allsamples.combined.filtered.vcf /mnt/projects/generic/data/ncbi/common_no_known_medical_impact_20140826.vcf
	PWD=$(pwd)
	(cd /data_synology/software/snpeff-3.6; java -jar SnpSift.jar annotate \
		-v /mnt/projects/generic/data/ncbi/common_no_known_medical_impact_20140826.vcf \
		<(cat /mnt/projects/p2ry8-crlf2/results/$< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@

pindel/allsamples.combined.filtered.dbsnp.snpeff.vcf: pindel/allsamples.combined.filtered.dbsnp.vcf
	PWD=$(pwd)
	(cd /data_synology/software/snpeff-3.6; java -Xmx4g -jar snpEff.jar -v -lof GRCh37.75 -stats $(PWD)/snpeff/$*.snpeff.summary.html $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

pindel/allsamples.combined.filtered.dbsnp.snpeff.dbNSFP.vcf: pindel/allsamples.combined.filtered.dbsnp.snpeff.vcf
	PWD=$(pwd)
	(cd /data_synology/software/snpeff-3.6; java -jar SnpSift.jar dbnsfp -v /mnt/projects/generic/data/dbNSFP-2.6/dbNSFP2.6_variant.tsv.gz -collapse \
		-f SIFT_pred,SIFT_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_score,SiPhy_29way_logOdds,LRT_pred,LRT_score,MutationTaster_pred,MutationTaster_score,MutationAssessor_pred,MutationAssessor_score,FATHMM_pred,FATHMM_score,RadialSVM_pred,RadialSVM_score,GERP++_RS,1000Gp1_AF,1000Gp1_AFR_AF,1000Gp1_EUR_AF,1000Gp1_AMR_AF,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF,Uniprot_acc,Interpro_domain, \
		$(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

pindel/allsamples.pindel.tsv: pindel/allsamples.combined.filtered.dbsnp.snpeff.dbNSFP.vcf /mnt/projects/p2ry8-crlf2/scripts/pindel-vcf-to-tsv.pl
	perl /mnt/projects/p2ry8-crlf2/scripts/pindel-vcf-to-tsv.pl \
		--vcf-in $< \
		--evs-file /mnt/projects/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		> $@.part
	mv $@.part $@ 

#---------------
# XENOME
#---------------
.PHONY: xenome
xenome: $(foreach P, $(PATIENTS_XENO), xenome/$P.xenome.txt)
xenome/%.xenome.txt: /mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_%_realigned.bam
	mkdir -p xenome
	(cd /home/cf/tools/xenome ; ./xenome classify \
		-v -T 8 -P idx \
		--dont-write-reads \
		-i <(java -jar /data_synology/software/picard-tools-1.114/SamToFastq.jar \
				VALIDATION_STRINGENCY=SILENT \
				INPUT=$<
				FASTQ=/dev/stdout \
				SECOND_END_FASTQ=/dev/null) \
		--output-filename-prefix $* \
	) 2>&1 1>$@.part | $(LOG)
	mv $@.part $@

#-----------	
# SNP PROFILE
#-----------	
	
.PHONY: snp
snp: snp-profile/allsamples.loh-segments.tsv

snp-profile/allsamples.loh-segments.tsv: $(foreach P, $(PATIENTS_MATCHED), snp-profile/$PD.snp-profile.pdf snp-profile/$PR.snp-profile.pdf snp-profile/$PC.snp-profile.pdf) \
	 									 $(foreach P, $(PATIENTS_DIA_ONLY), snp-profile/$PD.snp-profile.pdf snp-profile/$PC.snp-profile.pdf) \
	 									 $(foreach P, $(PATIENTS_REL2), snp-profile/$PR2.snp-profile.pdf snp-profile/$PC.snp-profile.pdf) \
	 									 $(foreach P, $(PATIENTS_REL3), snp-profile/$PR3.snp-profile.pdf snp-profile/$PC.snp-profile.pdf) \
	 									 $(foreach P, $(PATIENTS_XENO), snp-profile/$P.xeno.snp-profile.pdf)
	{ head -n1 snp-profile/839D.loh-segments.tsv; for f in snp-profile/*.loh-segments.tsv; do tail -n+2 "$$f"; done; } > $@.part
	mv $@.part $@
		
.SECONDEXPANSION:
snp-profile/%.snp-profile.pdf: coverage-genome/%.coverage-genome.tsv coverage-genome/$$(subst C.loh,C,$$(subst D.loh,C,$$(subst R.loh,C,$$(subst R3.loh,C,$$(subst R2.loh,C,%.loh))))).coverage-genome.tsv \
                               varscan/%.varscan.dbsnp.vcf varscan/$$(subst C.loh,C,$$(subst D.loh,C,$$(subst R.loh,C,$$(subst R3.loh,C,$$(subst R2.loh,C,%.loh))))).varscan.dbsnp.vcf \
                               /mnt/projects/generic/scripts/snp-profile.R
	mkdir -p snp-profile
	Rscript /mnt/projects/generic/scripts/snp-profile.R --sample-id $* --tumor-coverage $(word 1, $+) --normal-coverage $(word 2, $+) --tumor-vcf $(word 3, $+) --normal-vcf $(word 4, $+) --plot-output-file $@.part --loh-segments-output-file snp-profile/$*.loh-segments.tsv.part
	mv $@.part $@
	mv snp-profile/$*.loh-segments.tsv.part snp-profile/$*.loh-segments.tsv

snp-profile/%.xeno.snp-profile.pdf: coverage-genome/%.coverage-genome.tsv coverage-genome/948C.coverage-genome.tsv varscan/%.varscan.dbsnp.vcf /mnt/projects/generic/scripts/snp-profile.R
	mkdir -p snp-profile
	Rscript /mnt/projects/generic/scripts/snp-profile.R --sample-id $* --tumor-coverage $(word 1, $+) --normal-coverage $(word 2, $+) --tumor-vcf $(word 3, $+) --disomic-chromosome chr2 --plot-output-file $@.part --loh-segments-output-file snp-profile/$*.xeno.loh-segments.tsv.part
	mv $@.part $@
	mv snp-profile/$*.xeno.loh-segments.tsv.part snp-profile/$*.xeno.loh-segments.tsv

varscan/%.varscan.vcf: /mnt/projects/p2ry8-crlf2/data/bam/variant_calling_process_sample_%_realigned.bam
	mkdir -p varscan
	# -q 10 : remove non-uniquely mapping reads and reads with low mapping quality
	# -F 1024 : remove PCR and optical duplicates
	# -f 2 : only reads mapped in proper pairs
	java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar /data_synology/software/varscan-2.3.6/VarScan.v2.3.6.jar mpileup2cns \
		<(/data_synology/software/samtools-0.1.19/samtools view -u -q 10 -F 1024 -f 2 $(word 1,$^) | /data_synology/software/samtools-0.1.19/samtools mpileup -f /mnt/projects/generic/data/broad/human_g1k_v37.fasta -) \
		--min-coverage 8 \
		--min-reads2 4 \
		--min-avg-qual 20 \
		--min-var-freq 0.1 \
		--min-freq-for-hom 0.75 \
		--p-value 1 \
		--strand-filter 1 \
		--output-vcf 1 \
		--variants 1 \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@

varscan/%.varscan.dbsnp.vcf: varscan/%.varscan.vcf /mnt/projects/generic/data/ncbi/common_no_known_medical_impact_20140826.vcf
	mkdir -p snpeff
	PWD=$(pwd)
	(cd /data_synology/software/snpeff-3.6; java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar SnpSift.jar annotate \
		-v /mnt/projects/generic/data/ncbi/common_no_known_medical_impact_20140826.vcf \
		<(cat $(PWD)/$< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@
	rm $<

# SNP comparison (QC against sample swaps)
snp-comparison/allsamples.snp-comparison.tsv: $(foreach P, $(PATIENTS_MATCHED), snp-comparison/$PD.snp-comparison.tsv snp-comparison/$PR.snp-comparison.tsv) \
	 										  $(foreach P, $(PATIENTS_DIA_ONLY), snp-comparison/$PD.snp-comparison.tsv) \
	 										  $(foreach P, $(PATIENTS_REL2), snp-comparison/$PR2.snp-comparison.tsv) \
	 										  $(foreach P, $(PATIENTS_REL3), snp-comparison/$PR3.snp-comparison.tsv)
	rm -f $@
	{ head -n1 snp-comparison/839D.snp-comparison.tsv; for f in snp-comparison/*.snp-comparison.tsv; do tail -n+2 "$$f"; done; } > $@.part
	sort -k 3,3n $@.part > $@.sorted
	mv $@.sorted $@
	rm $@.part

.SECONDEXPANSION:
snp-comparison/%.snp-comparison.tsv: varscan/%.varscan.dbsnp.vcf varscan/$$(subst C.snp,C,$$(subst D.snp,C,$$(subst R.snp,C,$$(subst R3.snp,C,$$(subst R2.snp,C,%.snp))))).varscan.dbsnp.vcf /mnt/projects/generic/scripts/snp-comparison.R
	mkdir -p snp-comparison
	Rscript /mnt/projects/generic/scripts/snp-comparison.R --sample-id $* --tumor-vcf $(word 1, $+) --normal-vcf $(word 2, $+) --output-file $@.part
	mv $@.part $@
