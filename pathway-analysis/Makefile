export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'
DATA=/data/christian/p2ry8-crlf2/data/current

PATIENTS_TEST = 92
PATIENTS_MATCHED = 108 839 92 B36 BB16 GI13 HV57 HV80 LU3 N7 S23 SN18 737
PATIENTS_DIA_ONLY = 242 360 365 379 400 506 769 833 948
PATIENTS_REL2 = 108 737 
PATIENTS_REL3 = 715 

WIGFILES=$(foreach P, $(PATIENTS_MATCHED), $(DATA)/wig/$PD.wig $(DATA)/wig/$PC.wig $(DATA)/wig/$PR.wig) $(foreach P, $(PATIENTS_DIA_ONLY), $(DATA)/wig/$PD.wig $(DATA)/wig/$PC.wig) $(foreach P, $(PATIENTS_REL2), $(DATA)/wig/$PR2.wig) $(foreach P, $(PATIENTS_REL3), $(DATA)/wig/$PR3.wig)

all: pathways-pathwaycommons.tsv maf/relapsing.dia.maf dia-relapsing/smg.tsv rel/smg.tsv pathscan/sm_pathways.rel.annotated.tsv pathscan/sm_pathways.dia-relapsing.annotated.tsv

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

maf/relapsing.dia.maf: $(foreach P, $(PATIENTS_MATCHED), maf/$P_rem_dia.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1> $@.part
	cat $^ >> $@.part
	mv $@.part $@

maf/relapsing.rel.maf: $(foreach P, $(PATIENTS_MATCHED), maf/$P_rem_rel.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1> $@.part
	cat $^ >> $@.part
	mv $@.part $@

maf/%.maf: ../filtered_variants/%.dbsnp.snpeff.dbNSFP.filtered.vcf ../filtered_variants/%.indel.dbsnp.snpeff.dbNSFP.filtered.vcf ucsc-genes.g1k_v37.roi.gz ucsc-genes.g1k_v37.roi.gz.tbi ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv ~/p2ry8-crlf2/scripts/pathway-analysis/vcf2maf.pl
	mkdir -p maf
	cat $(word 1,$^) | perl ~/p2ry8-crlf2/scripts/pathway-analysis/vcf2maf.pl \
		--sample-tumor $(shell echo $* | perl -ne 's/_rem_dia/D/; s/_rem_rel/R/; print $$_;') \
		--sample-normal $(shell echo $* | perl -ne 's/_rem_dia/C/; s/_rem_rel/C/; print $$_;') \
		--music-roi ucsc-genes.g1k_v37.roi.gz \
		--mapping-entrez ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.20 \
		--deleterious \
		2>&1 1>$@.part | grep -vP '(Leading or trailing space|variant.Format|Domain annotations come from)'
	cat $(word 2,$^) | perl ~/p2ry8-crlf2/scripts/pathway-analysis/vcf2maf.pl \
		--sample-tumor $(shell echo $* | perl -ne 's/_rem_dia/D/; s/_rem_rel/R/; print $$_;') \
		--sample-normal $(shell echo $* | perl -ne 's/_rem_dia/C/; s/_rem_rel/C/; print $$_;') \
		--music-roi ucsc-genes.g1k_v37.roi.gz \
		--mapping-entrez ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.20 \
		--deleterious \
		2>&1 1>>$@.part | grep -vP '(Leading or trailing space|variant.Format|Domain annotations come from)'
	mv $@.part $@

$(DATA)/wig/%.wig: $(DATA)/bam/%.duplicate_marked.realigned.recalibrated.bam ~/generic/data/broad/human_g1k_v37.chromsizes.tsv 
	~/tools/bedtools-2.17.0/bin/genomeCoverageBed -bg -ibam $< -g ~/generic/data/broad/human_g1k_v37.chromsizes.tsv \
		| perl /home/STANNANET/christian.frech/hdall/scripts/bedgraph-to-wig.pl /dev/stdin $@.part
	mv $@.part $@

#--- pathway analysis diagnosis
	
dia-relapsing/total_covgs: ucsc-genes.g1k_v37.roi
	rm -f dia-relapsing/total_covgs
	rm -rf dia-relapsing/roi_covgs
	rm -rf dia-relapsing/gene_covgs
	rm -f dia-relapsing/wig-list.tsv
	for P in $(PATIENTS_MATCHED) ; do echo -e $$P"D\t$(DATA)/wig/"$$P"D.wig" >> dia-relapsing/wig-list.tsv ; done
	genome music bmr calc-wig-covg \
		--wig-list dia-relapsing/wig-list.tsv \
		--output-dir dia-relapsing/ \
		--reference-sequence ~/generic/data/broad/human_g1k_v37.fasta \
		--roi-file ucsc-genes.g1k_v37.roi

dia-relapsing/gene_mrs: dia-relapsing/total_covgs dia-relapsing/wig-list.tsv maf/relapsing.dia.maf ~/generic/data/broad/human_g1k_v37.fasta ucsc-genes.g1k_v37.roi
	rm -f dia-relapsing/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list dia-relapsing/wig-list.tsv \
		--maf-file maf/relapsing.dia.maf \
		--output-dir dia-relapsing/ \
		--reference-sequence ~/generic/data/broad/human_g1k_v37.fasta \
		--roi-file ucsc-genes.g1k_v37.roi

dia-relapsing/smg.tsv: dia-relapsing/gene_mrs
	genome music smg \
		--gene-mr-file dia-relapsing/gene_mrs \
		--output-file $@.part
	mv $@.part $@

pathscan/sm_pathways.dia-relapsing.raw.tsv: dia-relapsing/gene_mrs pathways-all.sizefiltered.tsv maf/relapsing.dia.maf rel/wig-list.tsv
	genome music path-scan \
		--bam-list dia-relapsing/wig-list.tsv \
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

#--- pathway analysis relapse

rel/total_covgs: ucsc-genes.g1k_v37.roi
	rm -f rel/total_covgs
	rm -rf rel/roi_covgs
	rm -rf rel/gene_covgs
	rm -f rel/wig-list.tsv
	for P in $(PATIENTS_MATCHED) ; do echo -e $$P"R\t$(DATA)/wig/"$$P"R.wig" >> rel/wig-list.tsv ; done
	genome music bmr calc-wig-covg \
		--wig-list rel/wig-list.tsv \
		--output-dir rel/ \
		--reference-sequence ~/generic/data/broad/human_g1k_v37.fasta \
		--roi-file ucsc-genes.g1k_v37.roi

rel/gene_mrs: rel/total_covgs rel/wig-list.tsv maf/relapsing.rel.maf ~/generic/data/broad/human_g1k_v37.fasta ucsc-genes.g1k_v37.roi
	rm -f rel/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list rel/wig-list.tsv \
		--maf-file maf/relapsing.rel.maf \
		--output-dir rel/ \
		--reference-sequence ~/generic/data/broad/human_g1k_v37.fasta \
		--roi-file ucsc-genes.g1k_v37.roi

rel/smg.tsv: rel/gene_mrs
	genome music smg \
		--gene-mr-file rel/gene_mrs \
		--output-file $@.part
	mv $@.part $@

pathscan/sm_pathways.rel.raw.tsv: rel/gene_mrs pathways-all.sizefiltered.tsv maf/relapsing.rel.maf rel/wig-list.tsv
	genome music path-scan \
		--bam-list rel/wig-list.tsv \
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
	R --no-save --quiet --slave -f ~/hdall/scripts/music/merge-sm-pathways.R \
		--args pathscan/sm_pathways.dia-relapsing.annotated.tsv pathscan/sm_pathways.rel.annotated.tsv \
		2>&1 1>$@.part | tee -a make.log
	R --no-save --quiet --slave -f ~/hdall/scripts/pathway-analysis/cluster-pathways.R \
		--args $@.part $@.clustered 1e-6 \
		2>&1 | tee -a make.log
	mv $@.clustered $@
	rm $@.part

