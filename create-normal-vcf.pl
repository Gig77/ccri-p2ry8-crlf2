use warnings FATAL => qw( all );
use strict;

use Vcf;

my $vcf = Vcf->new(file => "-");
$vcf->parse_header();
my (@samples) = $vcf->get_samples();

my $sample = $samples[0];
$sample = $samples[1] if ($sample !~ /C$/);
die "Not consitutional sample: $sample!\n" if ($sample !~ /C$/);

my ($patient) = $sample =~ /(.*)C$/;

while (my $x = $vcf->next_data_hash())
{
	next if ($x->{INFO}{SOMATIC});
		
	my $chr = $x->{CHROM};
	my $pos = $x->{POS};
	my $ref_allele = $x->{REF};
	my $alt_allele = $x->{ALT}[0];

    # INDEL VCF?
	if ($x->{INFO}{N_DP}) {
		##INFO=<ID=N_AC,Number=2,Type=Integer,Description="In NORMAL: # of reads supporting consensus indel/any indel at the site">
		##INFO=<ID=N_SC,Number=4,Type=Integer,Description="In NORMAL: strandness: counts of forward-/reverse-aligned indel-supporting reads / forward-/reverse-aligned reference supporting reads">
		##INFO=<ID=N_DP,Number=1,Type=Integer,Description="In NORMAL: total coverage at the site">
	
		#my ($n_sc_fwd_indel, $n_sc_rev_indel, $n_sc_fwd_ref, $n_sc_rev_ref) = split(",", $x->{INFO}{N_SC});
		#next if ($n_sc_fwd_indel == 0 or $n_sc_rev_indel == 0);
	
		my $n_dp = $x->{INFO}{N_DP};
		my ($n_ac) = split(",", $x->{INFO}{N_AC});
		
		next if ($n_ac < 2);
		
		print "$patient\t$chr\t$pos\t$ref_allele\t$alt_allele\t$n_dp\t$n_ac\tn/d\n";		
	}
	# SNP VCF
	else
	{
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		##FORMAT=<ID=FA,Number=A,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">
		my ($ad_ref, $ad_alt) = split(",", $x->{gtypes}{$sample}{AD});
		next if ($ad_alt < 2);

		my $gt = $x->{gtypes}{$sample}{GT};
		
		print "$patient\t$chr\t$pos\t$ref_allele\t$alt_allele\t".($ad_ref+$ad_alt)."\t$ad_alt\t$gt\n";		
	}
}
$vcf->close();
