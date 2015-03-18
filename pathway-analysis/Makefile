export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDEXPANSION:

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'
DATA=/data2/christian.frech/p2ry8-crlf2/data/current

PATIENTS_TEST = 92
PATIENTS_MATCHED = 108 92 715 737 839 B36 BB16 DL2 GI13 HV57 HV80 N7 S23 DS10898 SE15285 BJ17183 KE17247 AL9890 GL11356
PATIENTS_DIA_ONLY = 242 360 365 379 400 506 769 802 833 887 841 903 948 961 1060 1066 1089 HW11537 KT14158 TL14516
PATIENTS_REL2 = 108 737 
PATIENTS_REL3 = 715 737
PATIENTS_EXCLUDED_HYPER_REL = GI8 VS14645  # exclude hypermutated relapse samples
PATIENTS_EXCLUDED_PC_SUBCLONAL = MA5 BJ14367 LU3 SN18 460 545 564 957

all: pathways-pathwaycommons.tsv pathscan-dia-relapsing pathscan-rel pathscan/sm_pathways.dia-relapsing+rel.tsv

# --- prepare input pathways

~/generic/data/pathway-commons/PathwayCommons.4.All.GSEA.gmt: 
	curl http://www.pathwaycommons.org/pc2/downloads/Pathway%20Commons.4.All.GSEA.gmt.gz -o ~/generic/data/pathway-commons/PathwayCommons.4.All.GSEA.gmt.gz
	gunzip ~/generic/data/pathway-commons/PathwayCommons.4.All.GSEA.gmt.gz
	
pathways-pathwaycommons.tsv: ~/generic/data/pathway-commons/PathwayCommons.4.All.GSEA.gmt ~/p2ry8-crlf2/scripts/pathway-analysis/pathwaycommons2pathscan.pl 
	cat $< | perl ~/p2ry8-crlf2/scripts/pathway-analysis/pathwaycommons2pathscan.pl > $@.part
	mv $@.part $@ 

pathways-all.sizefiltered.tsv: pathways-pathwaycommons.tsv ~/hdall/results/music/pathscan/gene-sets.mutsigdb4.0.tsv
	cat $^ | perl -ne '@f = split /\t/; @g = split(/\|/, $$f[3]); print $$_ if (@g > 1 and @g < 400);' > $@.part
	mv $@.part $@


# --- prepare input mutations

maf/%.maf: ../filtered-variants/%.filtered.vcf ucsc-genes.g1k_v37.roi.gz ucsc-genes.g1k_v37.roi.gz.tbi ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv ~/p2ry8-crlf2/scripts/pathway-analysis/vcf2maf.pl
	mkdir -p maf
	cat $(word 1,$^) | perl ~/p2ry8-crlf2/scripts/pathway-analysis/vcf2maf.pl \
		--sample-tumor $(shell echo $* | perl -ne 's/_rem_dia/D/; s/_rem_rel/R/; print $$_;') \
		--sample-normal $(shell echo $* | perl -ne 's/_rem_dia/C/; s/_rem_rel/C/; print $$_;') \
		--music-roi ucsc-genes.g1k_v37.roi.gz \
		--mapping-entrez ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.20 \
		--deleterious \
		2>&1 1>$@.part | grep -vP '(Leading or trailing space|variant.Format|Domain annotations come from)'
	mv $@.part $@
	
#---------------------------------
#--- pathway analysis diagnosis, relapsing cases
#---------------------------------

.PHONY: pathscan-dia-relapsing
pathscan-dia-relapsing: dia-relapsing/smg.tsv pathscan/sm_pathways.dia-relapsing.annotated.tsv

maf/relapsing.dia.maf: $(foreach P, $(PATIENTS_MATCHED), maf/$P_rem_dia.maf)
	perl ~/p2ry8-crlf2/scripts/pathway-analysis/vcf2maf.pl --header 2>&1 1> $@.part
	cat $^ >> $@.part
	mv $@.part $@

dia-relapsing/bam-list.tsv: $(foreach P, $(PATIENTS_MATCHED), $(DATA)/bam/variant_calling_process_sample_$PD_realigned.bam $(DATA)/bam/variant_calling_process_sample_$PC_realigned.bam)
	mkdir -p dia-relapsing
	rm -f $@
	for P in $(PATIENTS_MATCHED) ; do echo -e $$P"D\t$(DATA)/bam/variant_calling_process_sample_"$$P"C_realigned.bam\t$(DATA)/bam/variant_calling_process_sample_"$$P"D_realigned.bam" >> $@.part ; done
	mv $@.part $@

dia-relapsing/roi_covgs/%.covg: $(DATA)/bam/variant_calling_process_sample_%_realigned.bam $(DATA)/bam/variant_calling_process_sample_$$(shell echo % | sed "s/[DR][23]\?$$$$//")C_realigned.bam ucsc-genes.g1k_v37.roi
	rm -rf dia-relapsing/tmp.$* ; mkdir -p dia-relapsing/tmp.$*
	echo -e "$*\t$(word 2, $^)\t$(word 1, $^)" > dia-relapsing/tmp.$*/bamlist
	genome music bmr calc-covg \
		--bam-list dia-relapsing/tmp.$*/bamlist \
		--output-dir dia-relapsing/tmp.$* \
		--reference-sequence ~/generic/data/broad/human_g1k_v37.fasta \
		--roi-file ucsc-genes.g1k_v37.roi
	mv dia-relapsing/tmp.$*/roi_covgs/* dia-relapsing/roi_covgs/
	mv dia-relapsing/tmp.$*/gene_covgs/* dia-relapsing/gene_covgs/
	rm -rf dia-relapsing/tmp.$*

dia-relapsing/total_covgs: $(foreach P, $(PATIENTS_MATCHED), dia-relapsing/roi_covgs/$PD.covg) dia-relapsing/bam-list.tsv ucsc-genes.g1k_v37.roi 
	rm -f dia-relapsing/total_covgs
	genome music bmr calc-covg \
		--bam-list dia-relapsing/bam-list.tsv \
		--output-dir dia-relapsing/ \
		--reference-sequence ~/generic/data/broad/human_g1k_v37.fasta \
		--roi-file ucsc-genes.g1k_v37.roi

dia-relapsing/gene_mrs: dia-relapsing/total_covgs dia-relapsing/bam-list.tsv maf/relapsing.dia.maf ~/generic/data/broad/human_g1k_v37.fasta ucsc-genes.g1k_v37.roi
	rm -f dia-relapsing/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list dia-relapsing/bam-list.tsv \
		--maf-file maf/relapsing.dia.maf \
		--output-dir dia-relapsing/ \
		--reference-sequence ~/generic/data/broad/human_g1k_v37.fasta \
		--roi-file ucsc-genes.g1k_v37.roi

dia-relapsing/smg.tsv: dia-relapsing/gene_mrs
	genome music smg \
		--gene-mr-file dia-relapsing/gene_mrs \
		--output-file $@.part
	mv $@.part $@

pathscan/sm_pathways.dia-relapsing.raw.tsv: dia-relapsing/gene_mrs pathways-all.sizefiltered.tsv maf/relapsing.dia.maf dia-relapsing/bam-list.tsv
	genome music path-scan \
		--bam-list dia-relapsing/bam-list.tsv \
		--gene-covg-dir dia-relapsing/gene_covgs/ \
		--maf-file maf/relapsing.dia.maf \
		--output-file $@.part \
		--pathway-file pathways-all.sizefiltered.tsv \
		--bmr $(shell grep Overall_BMR dia-relapsing/overall_bmrs | tail -1 | cut -f 4)
	mv $@.part $@
	mv $@.part_detailed $@_detailed

pathscan/sm_pathways.dia-relapsing.annotated.tsv: pathscan/sm_pathways.dia-relapsing.raw.tsv pathways-all.sizefiltered.tsv maf/relapsing.dia.maf ~/hdall/scripts/pathway-analysis/annotate-pathscan-result.pl
	perl ~/hdall/scripts/pathway-analysis/annotate-pathscan-result.pl \
		--pathway-file pathways-all.sizefiltered.tsv \
		--sm-pathways $< \
		--sm-pathways-detail $<_detailed \
		--maf-file maf/relapsing.dia.maf \
		> $@.part
	mv $@.part $@


#---------------------------------
#--- pathway analysis diagnosis, relapsing+nonrelapsing cases
# NOT IMPLEMENTED YET!
#---------------------------------

maf/dia.maf: $(foreach P, $(PATIENTS_MATCHED) $(PATIENTS_DIA_ONLY) $(PATIENTS_EXCLUDED_HYPER_REL), maf/$P_rem_dia.maf)
	perl ~/p2ry8-crlf2/scripts/pathway-analysis/vcf2maf.pl --header 2>&1 1> $@.part
	cat $^ >> $@.part
	mv $@.part $@

dia/bam-list.tsv: $(foreach P, $(PATIENTS_MATCHED) $(PATIENTS_DIA_ONLY) $(PATIENTS_EXCLUDED_HYPER_REL), $(DATA)/bam/variant_calling_process_sample_$PD_realigned.bam $(DATA)/bam/variant_calling_process_sample_$PC_realigned.bam)
	mkdir -p dia
	rm -f $@
	for P in $(PATIENTS_MATCHED) ; do echo -e $$P"D\t$(DATA)/bam/variant_calling_process_sample_"$$P"C_realigned.bam\t$(DATA)/bam/variant_calling_process_sample_"$$P"D_realigned.bam" >> $@.part ; done
	mv $@.part $@

#---------------------------------
#--- pathway analysis relapse
#---------------------------------

.PHONY: pathscan-rel
pathscan-rel: rel/smg.tsv pathscan/sm_pathways.rel.annotated.tsv 

maf/relapsing.rel.maf: $(foreach P, $(PATIENTS_MATCHED), maf/$P_rem_rel.maf)
	perl ~/p2ry8-crlf2/scripts/pathway-analysis/vcf2maf.pl --header 2>&1 1> $@.part
	cat $^ >> $@.part
	mv $@.part $@

rel/bam-list.tsv: $(foreach P, $(PATIENTS_MATCHED), $(DATA)/bam/variant_calling_process_sample_$PR_realigned.bam $(DATA)/bam/variant_calling_process_sample_$PC_realigned.bam)
	mkdir -p rel
	rm -f $@
	for P in $(PATIENTS_MATCHED) ; do echo -e $$P"R\t$(DATA)/bam/variant_calling_process_sample_"$$P"C_realigned.bam\t$(DATA)/bam/variant_calling_process_sample_"$$P"R_realigned.bam" >> $@.part ; done
	mv $@.part $@

rel/roi_covgs/%.covg: $(DATA)/bam/variant_calling_process_sample_%_realigned.bam $(DATA)/bam/variant_calling_process_sample_$$(shell echo % | sed "s/[DR][23]\?$$$$//")C_realigned.bam ucsc-genes.g1k_v37.roi
	rm -rf rel/tmp.$* ; mkdir -p rel/tmp.$*
	echo -e "$*\t$(word 2, $^)\t$(word 1, $^)" > rel/tmp.$*/bamlist
	genome music bmr calc-covg \
		--bam-list rel/tmp.$*/bamlist \
		--output-dir rel/tmp.$* \
		--reference-sequence ~/generic/data/broad/human_g1k_v37.fasta \
		--roi-file ucsc-genes.g1k_v37.roi
	mv rel/tmp.$*/roi_covgs/* rel/roi_covgs/
	mv rel/tmp.$*/gene_covgs/* rel/gene_covgs/
	rm -rf rel/tmp.$*

rel/total_covgs: $(foreach P, $(PATIENTS_MATCHED), rel/roi_covgs/$PR.covg) rel/bam-list.tsv ucsc-genes.g1k_v37.roi 
	rm -f rel/total_covgs
	genome music bmr calc-covg \
		--bam-list rel/bam-list.tsv \
		--output-dir rel/ \
		--reference-sequence ~/generic/data/broad/human_g1k_v37.fasta \
		--roi-file ucsc-genes.g1k_v37.roi

rel/gene_mrs: rel/total_covgs rel/bam-list.tsv maf/relapsing.rel.maf ~/generic/data/broad/human_g1k_v37.fasta ucsc-genes.g1k_v37.roi
	rm -f rel/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list rel/bam-list.tsv \
		--maf-file maf/relapsing.rel.maf \
		--output-dir rel/ \
		--reference-sequence ~/generic/data/broad/human_g1k_v37.fasta \
		--roi-file ucsc-genes.g1k_v37.roi

rel/smg.tsv: rel/gene_mrs
	genome music smg \
		--gene-mr-file rel/gene_mrs \
		--output-file $@.part
	mv $@.part $@

pathscan/sm_pathways.rel.raw.tsv: rel/gene_mrs pathways-all.sizefiltered.tsv maf/relapsing.rel.maf rel/bam-list.tsv
	genome music path-scan \
		--bam-list rel/bam-list.tsv \
		--gene-covg-dir rel/gene_covgs/ \
		--maf-file maf/relapsing.rel.maf \
		--output-file $@.part \
		--pathway-file pathways-all.sizefiltered.tsv \
		--bmr $(shell grep Overall_BMR rel/overall_bmrs | tail -1 | cut -f 4)
	mv $@.part $@
	mv $@.part_detailed $@_detailed

pathscan/sm_pathways.rel.annotated.tsv: pathscan/sm_pathways.rel.raw.tsv pathways-all.sizefiltered.tsv maf/relapsing.rel.maf ~/hdall/scripts/pathway-analysis/annotate-pathscan-result.pl
	perl ~/hdall/scripts/pathway-analysis/annotate-pathscan-result.pl \
		--pathway-file pathways-all.sizefiltered.tsv \
		--sm-pathways $< \
		--sm-pathways-detail $<_detailed \
		--maf-file maf/relapsing.rel.maf \
		> $@.part
	mv $@.part $@
	
#--- combine results into single list, cluster pathways by similarity

pathscan/sm_pathways.dia-relapsing+rel.tsv: pathscan/sm_pathways.dia-relapsing.annotated.tsv pathscan/sm_pathways.rel.annotated.tsv ~/hdall/scripts/pathway-analysis/cluster-pathways.R ~/hdall/scripts/music/merge-sm-pathways.R
	Rscript ~/hdall/scripts/music/merge-sm-pathways.R \
		pathscan/sm_pathways.dia-relapsing.annotated.tsv \
		pathscan/sm_pathways.rel.annotated.tsv \
		2>&1 1>$@.part | tee -a make.log
	Rscript ~/hdall/scripts/pathway-analysis/cluster-pathways.R $@.part $@.clustered 1e-6 \
		2>&1 | tee -a make.log
	mv $@.clustered $@
	rm $@.part

