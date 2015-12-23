use warnings FATAL => qw( all );
use strict;

use lib "/mnt/projects/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use List::Util qw(min max);
use Vcf;
use Data::Dumper;
use Getopt::Long;
use Tabix;
use Carp;

my ($vcf_in);
GetOptions
(
	"vcf-in=s" => \$vcf_in  # VCF input file
);

croak "ERROR: --vcf-in not specified" if (!$vcf_in);

my %rem_samples = (
	'108C' => 1,
	'92C' => 1,
	'460C' => 1,
	'545C' => 1,
	'564C' => 1,
	'715C' => 1,
	'737C' => 1,
	'839C' => 1,
	'B36C' => 1,
	'BB16C' => 1,
	'DL2C' => 1,
	'GI8C' => 1,
	'GI13C' => 1,
	'HV57C' => 1,
	'HV80C' => 1,
	'LU3C' => 1,
	'MA5C' => 1,
	'N7C' => 1,
	'S23C' => 1,
	'SN18C' => 1,
	'DS10898C' => 1,
	'VS14645C' => 1,
	'SE15285C' => 1,
	'BJ17183C' => 1,
	'KE17247C' => 1,
	'242C' => 1,
	'360C' => 1,
	'365C' => 1,
	'379C' => 1,
	'400C' => 1,
	'506C' => 1,
	'769C' => 1,
	'802C' => 1,
	'833C' => 1,
	'887C' => 1,
	'841C' => 1,
	'903C' => 1,
	'948C' => 1,
	'957C' => 1,
	'961C' => 1,
	'1060C' => 1,
	'1066C' => 1,
	'1089C' => 1,
	'HW11537C' => 1,
	'KT14158C' => 1,
	'TL14516C' => 1
);


$| = 1; # turn on autoflush

INFO("Processing file $vcf_in...");

my $vcf = Vcf->new(file => "$vcf_in");
$vcf->parse_header();

my (@samples) = $vcf->get_samples();

# copy header to stdout
print $vcf->format_header();

my %variants;
while (my $line = $vcf->next_line())
{
	my $x = $vcf->next_data_hash($line);

	my ($num_tum, $num_rem, $max_alt_rem, $max_alt_tum) = (0, 0, 0, 0);
	foreach my $s (@samples) 
	{
		next if ($s eq "715R3"); # 715 is actually HD ALL patient happened to be sequenced with Maria's cohort
		
		my $gt = $x->{gtypes}{$s}{GT};
		my ($ad_ref, $ad_alt) = split(",", $x->{gtypes}{$s}{AD});
		$ad_alt = $ad_ref if (!defined $ad_alt);

		if ($rem_samples{$s}) # remission sample?
		{
			$num_rem ++ if ($ad_alt >= 1);
			$max_alt_rem = max($max_alt_rem, $ad_alt); 
		}
		else
		{
			# check matched remission
			my ($rs) = $s =~ /(.*)(D|R\d?)$/;
			$rs .= "C";
			die "ERROR: Matched remission sample $rs not found\n" if (!defined $x->{gtypes}{$rs}{AD}); 
			my ($ad_ref_rem, $ad_alt_rem) = split(",", $x->{gtypes}{$rs}{AD});
			$ad_alt_rem = $ad_ref_rem if (!defined $ad_alt_rem);
			
			$num_tum ++ if ($ad_alt >= 6 and $ad_alt_rem == 0); # present in tumor but not matched remission
			$max_alt_tum = max($max_alt_tum, $ad_alt); 
		}
	}

	#print "$num_tum\t$num_rem\t$max_alt_tum\t$max_alt_rem\n" if ($num_tum > 0);
	
	next if ($num_tum == 0 or $num_rem > 1);

	print "$line";
}
	