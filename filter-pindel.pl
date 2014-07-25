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

my ($vcf_in);
GetOptions
(
	"vcf-in=s" => \$vcf_in  # VCF input file
);

croak "ERROR: --vcf-in not specified" if (!$vcf_in);

my %rem_samples = (
	'839C' => 1,
	'92C' => 1,
	'B36C' => 1,
	'BB16C' => 1,
	'GI13C' => 1,
	'HV57C' => 1,
	'HV80C' => 1,
	'LU3C' => 1,
	'N7C' => 1,
	'S23C' => 1,
	'SN18C' => 1,
	'242C' => 1,
	'360C' => 1,
	'365C' => 1,
	'379C' => 1,
	'400C' => 1,
	'506C' => 1,
	'769C' => 1,
	'833C' => 1,
	'948C' => 1,
	'737C' => 1,
	'108C' => 1
);


$| = 1; # turn on autoflush

INFO("Processing file $vcf_in...");

my $vcf = Vcf->new(file => "$vcf_in");
$vcf->parse_header();

my (@samples) = $vcf->get_samples();

# copy header
my $cmd = "grep -P '^#' $vcf_in";
system($cmd) == 0 or die "ERROR: grep vcf header failed: $cmd\n";

my %variants;
while (my $line = $vcf->next_line())
{
	my $x = $vcf->next_data_hash($line);

	my ($num_tum, $num_rem, $max_alt_rem, $max_alt_tum) = (0, 0, 0, 0);
	foreach my $s (@samples) 
	{
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
			$num_tum ++ if ($ad_alt >= 6);
			$max_alt_tum = max($max_alt_tum, $ad_alt); 
		}
	}

	next if ($num_rem > 0 or $num_tum < 1);

	print "$line";
}
	