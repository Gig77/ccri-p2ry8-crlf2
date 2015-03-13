use warnings FATAL => qw( all );
use strict;

use lib "$ENV{HOME}/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use List::Util qw(min max);
use Vcf;
use Data::Dumper;
use Getopt::Long;
use Tabix;
use Carp;

my ($vcf_out, $header, $rejected_variants_file, $sample_identifier, $vcf_in, $min_num_rem);
my ($rmsk_file, $simplerepeat_file, $blacklist_file, $segdup_file, $g1k_accessible_file, $ucsc_retro_file, $remission_variants_file, $evs_file);
GetOptions
(
	"sample=s" => \$sample_identifier, # e.g. 314_rem_dia
	"vcf-in=s" => \$vcf_in,  # VCF input file
	"vcf-out=s" => \$vcf_out,  # filtered VCF output file
	"header" => \$header,  # if set, write header line to output
	"rmsk-file=s" => \$rmsk_file, # TABIX indexed UCSC table rmsk
	"simpleRepeat-file=s" => \$simplerepeat_file, # TABIX indexed UCSC table rmsk
	"blacklist-file=s" => \$blacklist_file, # TABIX indexed UCSC table wgEncodeDacMapabilityConsensusExcludable
	"segdup-file=s" => \$segdup_file, # TABIX indexed UCSC table genomicSuperDups
	"g1k-accessible=s" => \$g1k_accessible_file, # TABIX indexed UCSC table tgpPhase1AccessibilityPilotCriteria
	"ucscRetro=s" => \$ucsc_retro_file, # TABIX indexed UCSC table ucscRetroAli5
	"rejected-variants-file=s" => \$rejected_variants_file, # file with variants rejected based on manual curation; will be filtered from output
	"remission-variants-file=s" => \$remission_variants_file, # TABIX indexed file with variants found in remission samples (GATK)
	"evs-file=s" => \$evs_file, # TABIX indexed file with wariants from Exome Variant Server (http://evs.gs.washington.edu/EVS/)
	"min-num-rem-to-exclude=i" => \$min_num_rem # minimum number of remission samples in which variant needs to be found in order to exclude it
);

$min_num_rem = 2 if (!defined $min_num_rem);

# TABLE: filtered-variants
if ($header)
{
	print "patient\t";		
	print "sample\t";
	print "cohort\t";
	print "var_type\t";
	print "status\t";
	print "rejected_because\t";
	print "chr\t";
	print "pos\t";
	print "dbSNP\t";
	print "ref\t";
	print "alt\t";
	print "gene\t";
	print "add_genes\t";
	print "impact\t";
	print "effect\t";
#	print "effect_snpeff\t";
	print "non_silent\t";
	print "deleterious\t";
	print "exons\t";
	print "dp_rem_tot\t";
	print "dp_rem_ref\t";
	print "dp_rem_var\t";
	print "freq_rem\t";
	print "dp_leu_tot\t";
	print "dp_leu_ref\t";
	print "dp_leu_var\t";
	print "freq_leu\t";
	print "aa_change\t";
	print "snpeff_effect\t";
	print "Polyphen2\t";
	print "SIFT\t";
	print "GERP++\t";
	print "SiPhy\t";
	print "InterPro\t";
	print "AF_1000G\t";
	print "repeat\t";
	print "segdup\t";
	print "blacklist\t";
	print "g1k-accessible\t";	
	print "retro\t";
	print "rem_samples\t";	
	print "evs_variant\t";
	print "AF_ESP6500\n";	
	exit;
}

my $debug = 1;

croak "ERROR: --sample not specified" if (!$sample_identifier);
croak "ERROR: --vcf-in not specified" if (!$vcf_in);
croak "ERROR: --rmsk-file not specified" if (!$rmsk_file);
croak "ERROR: --simpleRepeat-file not specified" if (!$simplerepeat_file);
croak "ERROR: --blacklist-file not specified" if (!$blacklist_file);
croak "ERROR: --segdup-file not specified" if (!$segdup_file);
croak "ERROR: --g1k-accessible not specified" if (!$g1k_accessible_file);
croak "ERROR: --ucscRetro not specified" if (!$ucsc_retro_file);
croak "ERROR: --remission-variants-file not specified" if (!$remission_variants_file);
croak "ERROR: --evs-file not specified" if (!$evs_file);

#-----------------------------------------------------------------------
# deduce VCF sample IDs from sample identifier provided by command line
#-----------------------------------------------------------------------

my %patient2sample = (
	'108_rem' => '108C',					'108_dia' => '108D',					'108_rel' => '108R1',			'108_rel2' => '108R2',
	'715_rem' => '715_Remission',			'715_dia' => '715_Diagnosis',			'715_rel' => '715_Relapse',		'715_rel3' => '715R3',						'm1957-715-rel_xeno' => 'M1957_715_Relapse1',
	'737_rem' => '737C',					'737_dia' => '737D',					'737_rel' => '737R',			'737_rel2' => '737R2',						'737_rel3' => '737_Relapse3',

	'92_rem' => '92C',						'92_dia' => '92D',						'92_rel' => '92R',
	'460_rem' => '460_Remission',			'460_dia' => '460_Diagnosis',			'460_rel' => '460_Relapse',
	'545_rem' => '545_Remission',			'545_dia' => '545_Diagnosis',			'545_rel' => '545_Relapse',		'm1963-545-rel_xeno' => 'm1963_545_Relapse',		'm1964-545-rel_xeno' => 'm1964_545Rn_Relapse',
	'564_rem' => '564_Remission',			'564_dia' => '564_Diagnosis',			'564_rel' => '564_Relapse',
	'839_rem' => '839C',					'839_dia' => '839D',					'839_rel' => '839R',
	'B36_rem' => 'B36C',					'B36_dia' => 'B36D',					'B36_rel' => 'B36R',
	'BB16_rem' => 'BB16C',					'BB16_dia' => 'BB16D',					'BB16_rel' => 'BB16R',
	'GI13_rem' => 'GI13C',					'GI13_dia' => 'GI13D',					'GI13_rel' => 'GI13R',
	'HV57_rem' => 'HV57C',					'HV57_dia' => 'HV57D',					'HV57_rel' => 'HV57R',
	'HV80_rem' => 'HV80C',					'HV80_dia' => 'HV80D',					'HV80_rel' => 'HV80R',
	'LU3_rem' => 'LU3C',					'LU3_dia' => 'LU3D',					'LU3_rel' => 'LU3R',
	'N7_rem' => 'N7C',						'N7_dia' => 'N7D',						'N7_rel' => 'N7R',
	'S23_rem' => 'S23C',					'S23_dia' => 'S23D',					'S23_rel' => 'S23R',
	'SN18_rem' => 'SN18C',					'SN18_dia' => 'SN18D',					'SN18_rel' => 'SN18R',
	'DL2_rem' => 'DL2_Remission',			'DL2_dia' => 'DL2_Diagnosis',			'DL2_rel' => '19981_DL2_R_Relapse',
	'GI8_rem' => 'GI8_Remission',			'GI8_dia' => 'GI8_Diagnosis',			'GI8_rel' => '19551_GI8_R_Relapse',
	'MA5_rem' => 'MA5_Remission',			'MA5_dia' => 'MA5_Diagnosis',			'MA5_rel' => '19319_MA5_R_Relapse',
	'VS14645_rem' => 'VS14645_Remission',	'VS14645_dia' => 'VS14645_Diagnosis',	'VS14645_rel' => '9931_Relapse',
	'BJ17183_rem' => 'BJ17183_Remission',	'BJ17183_dia' => 'BJ17183_Diagnosis',	'BJ17183_rel' => '14367_Relapse',
	'SE15285_rem' => 'SE1528_5_Remission',	'SE15285_dia' => 'SE15285_Diagnosis',	'SE15285_rel' => '13977_Relapse',
	'DS10898_rem' => 'DS10898_Remission',	'DS10898_dia' => 'DS10898_Diagnosis',	'DS10898_rel' => '5143_Relapse',
	'KE17247_rem' => 'KE17247_Remission',	'KE17247_dia' => 'KE17247_Diagnosis',	'KE17247_rel' => '15721_Relapse',

	'242_rem' => '242C',					'242_dia' => '242D',
	'360_rem' => '360C',					'360_dia' => '360D',
	'365_rem' => '365C',					'365_dia' => '365D',
	'379_rem' => '379C',					'379_dia' => '379D',
	'400_rem' => '400C',					'400_dia' => '400D',
	'506_rem' => '506C',					'506_dia' => '506D',
	'769_rem' => '769C',					'769_dia' => '769D',
	'833_rem' => '833C',					'833_dia' => '833D',
	'841_rem' => '842_Remission',			'841_dia' => '841_Diagnosis',
	'948_rem' => '948C',					'948_dia' => '948D',
	'802_rem' => '802_Remission',			'802_dia' => '802_Diagnosis',
	'887_rem' => '887_Remission',			'887_dia' => '887_Diagnosis',
	'903_rem' => '903_Remission',			'903_dia' => '903_Diagnosis',
	'957_rem' => '957_Remission',			'957_dia' => '957_Diagnosis',
	'961_rem' => '961_Remission',			'961_dia' => '961_Diagnosis',
	'1060_rem' => '1060_Remission',			'1060_dia' => '1060_Diagnosis',
	'1066_rem' => '1066_Remission',			'1066_dia' => '1066_Diagnosis',
	'1089_rem' => '1089_Remission',			'1089_dia' => '1089_Diagnosis',
	'HW11537_rem' => 'HW11537_Remission',	'HW11537_dia' => 'HW11537_Diagnosis',
	'KT14158_rem' => 'KT14158_Remission',	'KT14158_dia' => 'KT14158_Diagnosis',
	'TL14516_rem' => 'TL14516_Remission',	'TL14516_dia' => 'TL14516_Diagnosis',
	'AL9890_rem' => 'AL9890_Remission',		'AL9890_dia' => 'AL9890_Diagnosis',		'AL9890_rel' => 'AL9890_Relapse',
	'GL11356_rem' => 'GL11356_Remisson',	'GL11356_dia' => 'GL11356_Diagnosis',	'GL11356_rel' => 'GL11356_Relapse',
	
	'G_rem' => 'G_Remission',				'm1977-G-dia_xeno' => 'm1977_G_Dx_Diagnosis',
	'Y_rem' => 'Y3767_Remission',			'm1967-Y-rel_xeno' => 'm1967_Y_Relapse'
);

my %patient2cohort = (
	'108' => 'relapsing', 
	'92' => 'relapsing', 
	'460' => 'relapsing', 
	'545' => 'relapsing', 'm1963-545-rel_xeno' => 'relapsing', 'm1964-545-rel_xeno' => 'relapsing',
	'564' => 'relapsing', 
	'715' => 'relapsing', 'm1957-715-rel_xeno' => 'relapsing', 
	'737' => 'relapsing', 
	'839' => 'relapsing', 
	'B36' => 'relapsing', 
	'BB16' => 'relapsing', 
	'DL2' => 'relapsing', 
	'GI8' => 'relapsing', 
	'GI13' => 'relapsing', 
	'HV57' => 'relapsing',
	'HV80' => 'relapsing',
	'LU3' => 'relapsing',
	'MA5' => 'relapsing',
	'N7' => 'relapsing',
	'S23' => 'relapsing',
	'SN18' => 'relapsing',
	'DS10898' => 'relapsing',
	'VS14645' => 'relapsing',
	'SE15285' => 'relapsing',
	'BJ17183' => 'relapsing',
	'KE17247' => 'relapsing',
	'AL9890' => 'relapsing',
	'GL11356' => 'relapsing',
	'1060' => 'relapsing',
	
	'242' => 'non-relapsing',
	'360' => 'non-relapsing',
	'365' => 'non-relapsing',
	'379' => 'non-relapsing',
	'400' => 'non-relapsing',
	'506' => 'non-relapsing',
	'769' => 'non-relapsing',
	'802' => 'non-relapsing',
	'833' => 'non-relapsing',
	'887' => 'non-relapsing',
	'841' => 'non-relapsing',
	'903' => 'non-relapsing',
	'948' => 'non-relapsing',
	'957' => 'non-relapsing',
	'961' => 'non-relapsing',
	'1066' => 'non-relapsing',
	'1089' => 'non-relapsing',
	'HW11537' => 'non-relapsing',
	'KT14158' => 'non-relapsing',
	'TL14516' => 'non-relapsing',
	
	'm1977-G-dia' => 'HDALL',
	'm1967-Y-rel' => 'HDALL'
);

my ($patient, $vcf_sample_id_rem, $vcf_sample_id_tum) = $sample_identifier =~ /([^_]+)_([^_]+)_(.*)/ or croak "ERROR: could not parse sample identifier\n";
my $cmp_type = $vcf_sample_id_rem."_".$vcf_sample_id_tum;
die "ERROR: could not determine cohort for patient $patient\n" if (!$patient2cohort{$patient});
die "ERROR: invalid comparison type: $cmp_type\n" if ($cmp_type !~ /^(rem_dia|rem_rel\d?|rem_xeno)$/);

$vcf_sample_id_rem = $patient2sample{$patient."_$vcf_sample_id_rem"}; 
$vcf_sample_id_tum = $patient2sample{$patient."_$vcf_sample_id_tum"}; 

die "ERROR: Could not deduce normal VCF sample ID from provided sample identifier $sample_identifier\n" if (!$vcf_sample_id_rem);
die "ERROR: Could not deduce tumor VCF sample ID from provided sample identifier $sample_identifier\n" if (!$vcf_sample_id_tum);

print STDERR "Normal VCF sample ID: $vcf_sample_id_rem\n";
print STDERR "Tumor VCF sample ID: $vcf_sample_id_tum\n";

#-------------------------
# read in auxiliary files 
#-------------------------

my %rejected_variants;
if ($rejected_variants_file)
{
	open(R,"$rejected_variants_file") or die "could not open file $rejected_variants_file";
	<R>; # skip header
	while(<R>)
	{
		chomp;
		my ($patient, $sample, $var_type, $rejected_because, $chr, $pos) = split(/\t/);
		$rejected_variants{"$patient\t$sample\t$chr\t$pos"} = $rejected_because;
	}
	close(R);
	INFO(scalar(keys(%rejected_variants))." variants read from file $rejected_variants_file");
}

my $rmsk = Tabix->new(-data => $rmsk_file);
my $simpleRepeat = Tabix->new(-data => $simplerepeat_file);
my $blacklistdb = Tabix->new(-data => $blacklist_file);
my $segdup = Tabix->new(-data => $segdup_file);
my $g1kAcc = Tabix->new(-data => $g1k_accessible_file);
my $ucscRetro = Tabix->new(-data => $ucsc_retro_file);
my $remission = Tabix->new(-data => $remission_variants_file);
my $evs = Tabix->new(-data => $evs_file);

#-----------
# parse VCF 
#-----------

$| = 1; # turn on autoflush

my %variant_stati = 
(
	0 => 'wildtype',
	1 => 'germline',
	2 => 'somatic',
	3 => 'LOH',
	4 => 'post-transcriptional modification',
	5 => 'unknown'
);

INFO("Processing file $vcf_in...");

my $vcf = Vcf->new(file => "$vcf_in");
$vcf->parse_header();

my (@samples) = $vcf->get_samples();

# sanity checks
die "ERROR: Sample name $vcf_sample_id_rem not found!\n" if ($vcf_sample_id_rem ne $samples[0] and $vcf_sample_id_rem ne $samples[1]);
die "ERROR: Sample name $vcf_sample_id_tum not found!\n" if ($vcf_sample_id_tum ne $samples[0] and $vcf_sample_id_tum ne $samples[1]);
die "ERROR: Sample names identical: $vcf_sample_id_tum!\n" if ($vcf_sample_id_tum eq $vcf_sample_id_rem);

if ($vcf_out) 
{
	my $cmd = "grep -P '^#' $vcf_in > $vcf_out";
	system($cmd) == 0 or die "ERROR: grep vcf header failed: $cmd\n";
	open(VCFOUT,">>$vcf_out") or die "ERROR: could not write to file $vcf_out\n";
}


my ($tot_var, $filtered_qual, $filtered_gt, $filtered_alt, $filtered_germ) = (0, 0, 0, 0, 0);
my ($numrep, $num_blacklist, $numsegdup, $num_not_accessible, $num_retro, $num_remission, $num_evs) = (0, 0, 0, 0, 0, 0, 0);
my %qual_num;

while (my $line = $vcf->next_line())
{
	my $x = $vcf->next_data_hash($line);

	$tot_var ++;
	$qual_num{$x->{FILTER}->[0]} = $qual_num{$x->{FILTER}->[0]} ? $qual_num{$x->{FILTER}->[0]} + 1 : 1;
	
	if ($x->{gtypes}{$vcf_sample_id_rem}{GT} eq $x->{gtypes}{$vcf_sample_id_tum}{GT}) # no difference in genotype?
	{
		$filtered_gt ++;
		next;
	}
	
	my $gt_rem = $x->{gtypes}{$vcf_sample_id_rem}{GT};
	die "ERROR: Could not determine genotype of sample $vcf_sample_id_rem in file $vcf_in\n" if (!defined $gt_rem or $gt_rem eq "");

	if ($gt_rem =~ /1/) # germline variant?
	{
		$filtered_germ ++;
		next;
	}
	if (@{$x->{ALT}} != 1) # more than one alternative allele?
	{
		$filtered_alt ++;
		next;
	}		

	my $status = $x->{FILTER}->[0];
		
	if ($status eq "REJECT") # rejected by MuTect
	{
		$filtered_qual ++;
		next;
	}

	my ($dp_tum, $dp_rem, $freq_tum, $freq_rem, $ad_tum_ref, $ad_tum_alt, $ad_rem_ref, $ad_rem_alt);
	
#	print Dumper($x);
#	exit;		
	
	my $var_type;
	if (length($x->{REF}) == length($x->{ALT}->[0]))
	{
		$var_type = "snp";
		
		##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
		($ad_tum_ref, $ad_tum_alt) = split(",", $x->{gtypes}{$vcf_sample_id_tum}{AD});
		($ad_rem_ref, $ad_rem_alt) = split(",", $x->{gtypes}{$vcf_sample_id_rem}{AD});
		
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		($dp_tum, $dp_rem) = ($x->{gtypes}{$vcf_sample_id_tum}{DP}, $x->{gtypes}{$vcf_sample_id_rem}{DP});

		##FORMAT=<ID=FA,Number=A,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">
		$freq_tum = $x->{gtypes}{$vcf_sample_id_tum}{FA};
		$freq_rem = $x->{gtypes}{$vcf_sample_id_rem}{FA};
		
		#next if ($status eq "REJECT" and ($dp_tum <= 50 or $dp_rem <= 50 or $freq_tum < 0.2 or $freq_rem > 0.05));	
	}
	else
	{
		$var_type = "indel";
		
		##INFO=<ID=T_DP,Number=1,Type=Integer,Description="In TUMOR: total coverage at the site">
		##INFO=<ID=N_DP,Number=1,Type=Integer,Description="In NORMAL: total coverage at the site">
		($dp_tum, $dp_rem) = ($x->{INFO}{T_DP}, $x->{INFO}{N_DP});
		
		##INFO=<ID=T_AC,Number=2,Type=Integer,Description="In TUMOR: # of reads supporting consensus indel/any indel at the site">
		##INFO=<ID=N_AC,Number=2,Type=Integer,Description="In NORMAL: # of reads supporting consensus indel/any indel at the site">
		my ($ad_tum_any_indel, $ad_rem_any_indel);
		($ad_tum_alt, $ad_tum_any_indel) = split(",", $x->{INFO}{T_AC}); 		
		($ad_rem_alt, $ad_rem_any_indel) = split(",", $x->{INFO}{N_AC});
		($ad_tum_ref, $ad_rem_ref) = ($dp_tum - $ad_tum_any_indel, $dp_rem - $ad_rem_any_indel); 		
		$freq_tum = sprintf("%.3f", $ad_tum_alt / $dp_tum);
		$freq_rem = sprintf("%.3f", $ad_rem_alt / $dp_rem);

		# insufficient read depth
		if ($dp_tum < 10)
		{
			#INFO("REJECT: READ DEPTH < 10: $dp_tum");
			next;			
		}

		# require high consensus call for indel
		if ($ad_tum_alt/$ad_tum_any_indel < 0.7)
		{
			#INFO("REJECT: BAD CONSENSUS: ",$x->{CHROM},":",$x->{POS},"");
			next;
		}		
		
		# require support from both strands
		##INFO=<ID=T_SC,Number=4,Type=Integer,Description="In TUMOR: strandness: counts of forward-/reverse-aligned indel-supporting reads / forward-/reverse-aligned reference supporting reads">
		my ($reads_indel_fwd, $reads_indel_rev, $reads_ref_fwd, $reads_ref_rev) = split(",", $x->{INFO}{T_SC});
		if ($reads_indel_fwd == 0 or $reads_indel_rev == 0)
		{
			#INFO("REJECT: STRAND BIAS: ",$x->{CHROM},":",$x->{POS},"\t","reads_indel_fwd: $reads_indel_fwd\treads_indel_rev: $reads_indel_rev");
			next;
		}
		
		# check alignment quality around indel
		##INFO=<ID=T_NQSMM,Number=2,Type=Float,Description="In TUMOR: Within NQS window: fraction of mismatching bases in consensus indel-supporting reads/in reference-supporting reads">
		my ($frac_mm_reads_indel, $frac_mm_reads_ref) = split(",", $x->{INFO}{T_NQSMM});
		if ($frac_mm_reads_indel - $frac_mm_reads_ref > 0.01)
		{
			#INFO("REJECT: POOR ALIGNMENT: ",$x->{CHROM},":",$x->{POS},"\t","frac_mm_reads_indel: $frac_mm_reads_indel\tfrac_mm_reads_ref: $frac_mm_reads_ref");
			next;
		}

		# check mapping quality
		##INFO=<ID=T_MQ,Number=2,Type=Float,Description="In TUMOR: average mapping quality of consensus indel-supporting reads/reference-supporting reads">
		##INFO=<ID=N_MQ,Number=2,Type=Float,Description="In NORMAL: average mapping quality of consensus indel-supporting reads/reference-supporting reads">
		my ($mq_indel_tum, $mq_ref_tum) = split(",", $x->{INFO}{T_MQ});
		my ($mq_indel_rem, $mq_ref_rem) = split(",", $x->{INFO}{N_MQ});		
		if ($mq_indel_tum < 40)
		{
			#INFO("REJECT: POOR MAPPING: ",$x->{CHROM},":",$x->{POS},"\t","T_MQ=$mq_indel_tum,$mq_ref_tum");
			next;
		}
		
				
#		print "reads_indel_fwd: $reads_indel_fwd\n";
#		print "reads_indel_rev: $reads_indel_rev\n";
#		print "reads_ref_fwd: $reads_ref_fwd\n";
#		print "reads_ref_rev: $reads_ref_rev\n";
	}

	my (@repeats, @dups, @blacklist, @retro, @rem_samples, %evss);
	my ($chr, $pos) = ("chr".$x->{CHROM}, $x->{POS});

	# ----- annotate overlapping repeat regions
	{
		my $iter = $rmsk->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $rmsk->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@repeats, "$s[10]:$s[11]:$s[12]");
			}
		}		
	}

	{
		my $iter = $simpleRepeat->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $simpleRepeat->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@repeats, "$s[16]($s[6])");
			}
		}		
	}
	$numrep ++ if (@repeats > 0);

	# ----- annotate segmental duplications
	{
		my $iter = $segdup->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $segdup->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@dups, "$s[4]:".sprintf("%.2f", $s[26]));
			}		
		}		
		$numsegdup ++ if (@dups > 0);
	}	
	
	# ----- annotate overlapping DAC blacklisted regions
	{
		my $iter = $blacklistdb->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $blacklistdb->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@blacklist, "$s[4]");
			}		
		}
		$num_blacklist ++ if (@blacklist > 0);		
	}

	# ----- annotate overlapping g1k accessible regions
	my $accessible = "no";
	{
		my $iter = $g1kAcc->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $g1kAcc->read($iter)) 
			{
				$accessible = "";
				last;
			}		
		}
		$num_not_accessible ++ if ($accessible eq "no");		
	}

	# ----- annotate overlapping retrotransposed (pseudo) genes
	{
		my $iter = $ucscRetro->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $ucscRetro->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@retro, $s[10]);
			}		
		}
		$num_retro ++ if (@retro > 0);
	}

	# ----- annotate variants found in remission samples
	{
		my $iter = $remission->query($x->{CHROM}, $pos-1, $pos);
		if ($iter and $iter->{_})
		{
			while (my $line = $remission->read($iter)) 
			{
				my ($sample, $rchr, $rpos, $ref_allele, $alt_allele, $dp, $ad, $gt) = split("\t", $line);
				if ($pos eq $rpos and $x->{REF} eq $ref_allele and $x->{ALT}->[0] eq $alt_allele and $ad >= 3 and $ad/$dp > 0.05)
				{
					push(@rem_samples, "$sample($ad)");
				}
			}		
		}
		$num_remission ++ if (@rem_samples > 0);
	}

	# ----- annotate variants found in Exome Variant Server
	{
		my $iter = $evs->query($chr, $pos-1, $pos);
		if ($iter and $iter->{_})
		{
			while (my $line = $evs->read($iter)) 
			{
				my ($echr, $epos, $rsID, $dbSNPVersion, $alleles, $europeanAmericanAlleleCount, $africanAmericanAlleleCount, $allAlleleCount, $MAFinPercent_EA_AA_All, $europeanAmericanGenotypeCount, 
					$africanAmericanGenotypeCount, $allGenotypeCount, $avgSampleReadDepth, $genes, $geneAccession, $functionGVS, $hgvsProteinVariant, $hgvsCdnaVariant, $codingDnaSize, 
					$conservationScorePhastCons, $conservationScoreGERP, $granthamScore, $polyphen2_score, $refBaseNCBI37, $chimpAllele, $clinicalInfo, $filterStatus, $onIlluminaHumanExomeChip,
					$gwasPubMedInfo, $EA_EstimatedAge_kyrs, $AA_EstimatedAge_kyrs) = split(/\s/, $line);
					
				next if ($echr ne $chr or $epos ne $pos);
				foreach my $allele (split(";", $alleles))
				{
					my ($ref, $alt) = $allele =~ /(.+)\>(.+)/;
					if ($ref eq $x->{REF} and $alt eq $x->{ALT}->[0])
					{
						my ($alt_count, $ref_count) = $allAlleleCount =~ /(\d+).+?(\d+)/;
						my $alt_percent = sprintf("%.3f", $alt_count/($alt_count+$ref_count)*100);
						$evss{$alt_percent} = 1;					
					}
				}
			}		
		}
		$num_evs ++ if (keys(%evss) > 0);
	}
	
	# ----- G1K frequency
	my $g1k_max_freq = 0;
	foreach my $f ($x->{INFO}{'dbNSFP_1000Gp1_AF'}, $x->{INFO}{'dbNSFP_1000Gp1_AFR_AF'}, $x->{INFO}{'dbNSFP_1000Gp1_EUR_AF'}, $x->{INFO}{'dbNSFP_1000Gp1_AMR_AF'}, $x->{INFO}{'dbNSFP_1000Gp1_ASN_AF'})
	{
		$g1k_max_freq = $f if (defined $f and $g1k_max_freq < $f);
	}
	
	# ----- ESP6500 frequency
	my $evs_freq = (defined $x->{INFO}{dbNSFP_ESP6500_AA_AF} and defined $x->{INFO}{dbNSFP_ESP6500_EA_AF}) 
					? max($x->{INFO}{dbNSFP_ESP6500_AA_AF}, $x->{INFO}{dbNSFP_ESP6500_EA_AF}) 
					: defined $x->{INFO}{dbNSFP_ESP6500_AA_AF} 
						? $x->{INFO}{dbNSFP_ESP6500_AA_AF} 
						: defined $x->{INFO}{dbNSFP_ESP6500_EA_AF} 
							? $x->{INFO}{dbNSFP_ESP6500_EA_AF} 
							: 0;
	$num_evs ++ if ($evs_freq > 0);
	
	my @rejected_because;
	if ($rejected_variants{"$patient\t$cmp_type\t".$x->{CHROM}."\t".$x->{POS}}) { push(@rejected_because, "manual inspection (".$rejected_variants{"$patient\t$cmp_type\t".$x->{CHROM}."\t".$x->{POS}}.")")}
	if (@repeats > 0) { push(@rejected_because, "repetitive region"); }
	if (@dups > 0) { push(@rejected_because, "segmental duplication"); }
	if (@blacklist > 0) { push(@rejected_because, "blacklisted region"); }
	#if (@retro > 0) { push(@rejected_because, "retrotransposon"); }
	if (@rem_samples >= $min_num_rem) { push(@rejected_because, "present remissions"); }
	if ($x->{ID} and $x->{ID} ne ".")  { push(@rejected_because, "dbSNP"); }  
	if ($g1k_max_freq > 0.01) { push(@rejected_because, "G1K"); }
	if ($evs_freq > 0.01) { push(@rejected_because, "ESP6500"); }
	if ($x->{CHROM} eq "hs37d5") { push(@rejected_because, "decoy genome"); }
	
	# keep known pathogenic CRLF2 mutations overlapping with segmental duplication (PAR1 on Y chromosome!)
	my ($gene, $add_genes, $impact, $effect, $aa_change) = get_impact($x->{INFO}{EFF});
	@rejected_because = () if ($gene eq "CRLF2" && $aa_change eq "F232C"); 
	
	
	$line =~ s/^([^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t)[^\t]+/$1REJECT/ if (@rejected_because > 0);
	$line =~ s/^([^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t)[^\t]+/$1MISSED/ if ($status eq "MISSED");
	
	my $polyphen = $x->{INFO}{'dbNSFP_Polyphen2_HVAR_pred'};
	my $sift = undef;
	if (defined $x->{INFO}{'dbNSFP_SIFT_score'})
	{
		foreach my $s (split(",", $x->{INFO}{'dbNSFP_SIFT_score'}))
		{
			next if (!defined $s or $s eq ".");
			$sift = $s if (!defined $sift or $s < $sift);	
		}
	} 
	my $siphy = $x->{INFO}{'dbNSFP_SiPhy_29way_logOdds'};
	
	my $is_deleterious = "n/d";
	$is_deleterious = "yes" if ($effect eq "NON_SYNONYMOUS_CODING" and $polyphen and $polyphen =~ /D/ and defined $sift and $sift < 0.05); # polyphen & sift
	$is_deleterious = "yes" if ($effect eq "NON_SYNONYMOUS_CODING" and $polyphen and $polyphen =~ /D/ and defined $siphy and $siphy >= 12); # polyphen & siphy
	$is_deleterious = "yes" if ($effect eq "NON_SYNONYMOUS_CODING" and defined $sift and $sift < 0.05 and defined $siphy and $siphy >= 12); # sift and siphy
	$is_deleterious = "yes" if ($effect eq "NON_SYNONYMOUS_CODING" and defined $siphy and $siphy > 20); # siphy only, highly conserved (keeps GNAQ)
	$is_deleterious = "yes" if ($effect eq "FRAME_SHIFT" or $effect eq "SPLICE_SITE_ACCEPTOR" or $effect eq "SPLICE_SITE_DONOR" or $effect eq "STOP_GAINED");
	$is_deleterious = "no" if ($is_deleterious ne "yes" and defined $polyphen and defined $sift);
	$is_deleterious = "no" if ($effect eq "DOWNSTREAM" or $effect eq "UPSTREAM" or $effect eq "INTRON" or $effect eq "INTERGENIC" or $effect eq "SYNONYMOUS_CODING" or $effect eq "SYNONYMOUS_STOP" or $effect eq "SYNONYMOUS_START" or $effect eq "UTR_3_PRIME" or $effect eq "UTR_5_PRIME" or $effect eq "UTR_5_DELETED" or $effect eq "UTR_3_DELETED" or $effect eq "START_GAINED");
	
	my $non_silent = 0;
	$non_silent = 1 if ($effect eq "STOP_GAINED" or $effect eq "STOP_LOST" or $effect eq "SPLICE_SITE_DONOR" or $effect eq "SPLICE_SITE_ACCEPTOR" or $effect eq "FRAME_SHIFT" or $effect eq "CODON_CHANGE_PLUS_CODON_INSERTION" or $effect eq "CODON_DELETION" or $effect eq "NON_SYNONYMOUS_CODING" or $effect eq "CODON_INSERTION" or $effect eq "CODON_CHANGE_PLUS_CODON_DELETION" or $effect eq "NON_SYNONYMOUS_START" or $effect eq "START_LOST");
	
	print VCFOUT "$line" if ($vcf_out);
	
	print "$patient\t";		
	print "$cmp_type\t";
	print $patient2cohort{$patient}, "\t";
	print "$var_type\t";
	print @rejected_because > 0 ? "REJECT\t" : "$status\t";
	print @rejected_because > 0 ? join(";", @rejected_because) : "", "\t";
	print $x->{CHROM},"\t";
	print $x->{POS},"\t";
	print $x->{ID},"\t";
	print $x->{REF},"\t";
	print $x->{ALT}->[0],"\t";
	print "$gene\t";
	print "$add_genes\t";
	print "$impact\t";
#	print $x->{INFO}{SNPEFF_FUNCTIONAL_CLASS} ? $x->{INFO}{SNPEFF_FUNCTIONAL_CLASS} : "", "\t"; 
#	print $x->{INFO}{SNPEFF_IMPACT} ? $x->{INFO}{SNPEFF_IMPACT} : "", "\t"; 
	print "$effect\t";
#	print $x->{INFO}{SNPEFF_EFFECT} ? $x->{INFO}{SNPEFF_EFFECT} : "", "\t"; 
	print "$non_silent\t";
	print "$is_deleterious\t";
	print $x->{INFO}{SNPEFF_EXON_ID} ? $x->{INFO}{SNPEFF_EXON_ID} : "", "\t"; 
#	print join(",", @{$x->{FILTER}}),"\t";
#	print exists $x->{gtypes}{$vcf_sample_id_tum}{SS} ? $variant_stati{$x->{gtypes}{$vcf_sample_id_tum}{SS}} : "n/a", "\t";
	print "$dp_rem\t";
	print "$ad_rem_ref\t";
	print "$ad_rem_alt\t";
	print "$freq_rem\t";
	print "$dp_tum\t";
	print "$ad_tum_ref\t";
	print "$ad_tum_alt\t";
	print "$freq_tum\t";
	print "$aa_change\t";
	print "EFF=",$x->{INFO}{EFF},"\t";
	print defined $polyphen ? $polyphen : "", "\t"; # Polyphen2 prediction based on HumVar, 'D' ('porobably damaging'), 'P' ('possibly damaging') and 'B' ('benign'). Multiple entries separated by ';' 
	print defined $sift ? $sift : "", "\t"; # SIFT score, If a score is smaller than 0.05 the corresponding NS is predicted as 'D(amaging)'; otherwise it is predicted as 'T(olerated)'
	print defined $x->{INFO}{'dbNSFP_GERP++_RS'} ? $x->{INFO}{'dbNSFP_GERP++_RS'} : "", "\t"; # GERP++ RS score, the larger the score, the more conserved the site 
	print defined $siphy ? $siphy : "", "\t"; # SiPhy score based on 29 mammals genomes. The larger the score, the more conserved the site.
	my $domains = $x->{INFO}{'dbNSFP_Interpro_domain'}; # domain or conserved site on which the variant locates. Domain annotations come from Interpro database. The number in the brackets following a specific domain is the count of times Interpro assigns the variant position to that domain, typically coming from different predicting databases
	if ($domains)
	{
		$domains =~ s/\),/\)\|/g;
		$domains =~ s/\|$//;
		print "$domains\t";
	}
	else
	{
		print "\t";
	}
	print $g1k_max_freq > 0 ? $g1k_max_freq : "";  # alternative allele frequency in the whole 1000Gp1 data
	print "\t", join(',', @repeats), "\t", join(',', @dups), "\t", join(',', @blacklist), "\t$accessible\t", join(",", @retro), "\t", join(",", @rem_samples), "\t", join(";", keys(%evss));
	print "\t", $evs_freq > 0 ? $evs_freq : ""; 
	print "\n";
			
#	print "\n"; print Dumper($x); exit;
}

# add 737 KRAS conserved clonal variant missed by MuTect b/c three reads were present in remission sample, presumably because of MRD
if ($patient eq "737" and $cmp_type eq "rem_dia") {
	print("737\trem_dia\trelapsing\tsnp\tMISSED\t\t12\t25398284\t.\tC\tT\tKRAS\t\tMODERATE\tNON_SYNONYMOUS_CODING\t1\tyes\t2\t150\t147\t3\t0.02\t187\t115\t70\t0.374\tG12D\tEFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|188|KRAS|protein_coding|CODING|ENST00000311936|2|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|189|KRAS|protein_coding|CODING|ENST00000256078|2|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|43|KRAS|protein_coding|CODING|ENST00000556131|2|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|75|KRAS|protein_coding|CODING|ENST00000557334|2|1)\tB,B\t0\t5.68\t18.3719\tSmall_GTP-binding_protein_domain_(1)\t\t\t\t\t\t\t\t\t\n");
} elsif ($patient eq "737" and $cmp_type eq "rem_rel") {
	print("737\trem_rel\trelapsing\tsnp\tMISSED\t\t12\t25398284\t.\tC\tT\tKRAS\t\tMODERATE\tNON_SYNONYMOUS_CODING\t1\tyes\t2\t150\t147\t3\t0.02\t154\t81\t73\t0.474\tG12D\tEFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|188|KRAS|protein_coding|CODING|ENST00000311936|2|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|189|KRAS|protein_coding|CODING|ENST00000256078|2|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|43|KRAS|protein_coding|CODING|ENST00000556131|2|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|75|KRAS|protein_coding|CODING|ENST00000557334|2|1)\tB,B\t0\t5.68\t18.3719\tSmall_GTP-binding_protein_domain_(1)\t\t\t\t\t\t\t\t\t\n");		
} elsif ($patient eq "737" and $cmp_type eq "rem_rel2") {
	print("737\trem_rel2\trelapsing\tsnp\tMISSED\t\t12\t25398284\t.\tC\tT\tKRAS\t\tMODERATE\tNON_SYNONYMOUS_CODING\t1\tyes\t2\t150\t147\t3\t0.02\t152\t138\t14\t0.0921\tG12D\tEFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|188|KRAS|protein_coding|CODING|ENST00000311936|2|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|189|KRAS|protein_coding|CODING|ENST00000256078|2|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|43|KRAS|protein_coding|CODING|ENST00000556131|2|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|75|KRAS|protein_coding|CODING|ENST00000557334|2|1)\tB,B\t0\t5.68\t18.3719\tSmall_GTP-binding_protein_domain_(1)\t\t\t\t\t\t\t\t\t\n");		
} elsif ($patient eq "737" and $cmp_type eq "rem_rel3") {
	print("737\trem_rel3\trelapsing\tsnp\tMISSED\t\t12\t25398284\t.\tC\tT\tKRAS\t\tMODERATE\tNON_SYNONYMOUS_CODING\t1\tyes\t2\t150\t147\t3\t0.02\t141\t89\t52\t0.368794326\tG12D\tEFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|188|KRAS|protein_coding|CODING|ENST00000311936|2|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|189|KRAS|protein_coding|CODING|ENST00000256078|2|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|43|KRAS|protein_coding|CODING|ENST00000556131|2|1),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gGt/gAt|G12D|75|KRAS|protein_coding|CODING|ENST00000557334|2|1)\tB,B\t0\t5.68\t18.3719\tSmall_GTP-binding_protein_domain_(1)\t\t\t\t\t\t\t\t\t\n");		
}	

$vcf->close();
close(VCFOUT) if ($vcf_out);
	
if ($debug)
{
	INFO("  Total number of variants: $tot_var");
	INFO("  Variants by quality:");
	foreach my $k (keys(%qual_num))
	{
		INFO("    $k: ", $qual_num{$k});
	}
	INFO("  Rejected by MuTect: $filtered_qual");
	INFO("  Excluded due to equal genotype: $filtered_gt");
	INFO("  Excluded due to missing alternative allele: $filtered_alt");
	INFO("  Excluded germline variants: $filtered_germ");
	INFO("  $numrep variants annotated with overlapping repeat.");
	INFO("  $num_blacklist variants annotated with overlapping blacklisted region.");
	INFO("  $numsegdup variants annotated with overlapping segmental duplication.");
	INFO("  $num_not_accessible variants fall into G1K non-accessible region.");
	INFO("  $num_retro variants annotated with overlapping retrotransposed (pseudo)gene.");
	INFO("  $num_remission variants present in remission sample(s).");
	INFO("  $num_evs variants present in Exome Variant Server (EVS).");
}

# ------------------------------------------

sub get_impact
{
	my $effs = shift or die "ERROR: effect not specified";

	# determine all genes impacted by variants
	my (%genes_by_impact, %all_genes, $combined_effect, $combined_impact, %aa_changes);
	foreach my $eff (split(",", $effs))
	{
		my ($effect, $rest) = $eff =~ /([^\(]+)\(([^\)]+)\)/
			or die "ERROR: could not parse SNP effect: $effs\n";

		my ($impact, $class, $codon_change, $aa_change, $aa_length, $gene_name, $gene_biotype, 
			$coding, $transcript, $exon, $genotype_num) = split('\|', $rest)
				or die "ERROR: could not parse SNP effect: $eff\n";
		 
		$aa_changes{$aa_change} = 1 if ($aa_change);

		if ($exon and $transcript and $gene_name)
		{
			$transcript =~ s/\.\d+$//; # remove version number from accession
			$transcript =~ s/\.\d+$//; 
		}
			
		# gene impacted by variant?
		if ($gene_name)
		{
			$genes_by_impact{$impact}{$gene_name} = $effect;
			$all_genes{$gene_name} = 1;
		}

		$combined_impact = $impact;		
		$combined_effect = $effect;
	}
	
	# if multiple genes are affected, preferentially chose gene with the predicted higher impact
	if ($genes_by_impact{HIGH})
	{
		$combined_impact = "HIGH";
	}
	elsif ($genes_by_impact{MODERATE})
	{
		$combined_impact = "MODERATE";
	}
	elsif ($genes_by_impact{LOW})
	{
		$combined_impact = "LOW";
	}
	elsif ($genes_by_impact{MODIFIER})
	{
		$combined_impact = "MODIFIER";
	}
	
	my ($gene, $add_genes) = ("", "");
	if (keys(%all_genes) > 0)
	{
		my @sorted_genes = sort keys(%{$genes_by_impact{$combined_impact}});
		$gene = $sorted_genes[0]; # first choice is first in alphabetically sorted list
		if ($gene =~ /^LOC/) # if this is a generic gene name, try to find non-generic one instead
		{
			foreach my $g (@sorted_genes)
			{
				if ($g !~ /^LOC/)
				{
					$gene = $g;
					last;
				}	
			}
		}
		$combined_effect = $genes_by_impact{$combined_impact}{$gene};
		delete $all_genes{$gene};
		$add_genes = join(",", keys(%all_genes));
	}

	return ($gene, $add_genes, $combined_impact, $combined_effect, join(";", keys(%aa_changes)));
}
