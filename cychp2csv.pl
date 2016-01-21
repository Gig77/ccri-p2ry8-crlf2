use warnings FATAL => qw( all );
use strict;

my $region = $ARGV[0];
my ($chr, $start, $end);
if (defined $region) {
	($chr, $start, $end) = $region =~ /([^:]+):?(\d+)?\-?(\d+)?/;
	die "ERROR: Could not parse specified region: $region\n" if (!$chr && !$start && !$end);
	print STDERR "Parsing region $region...\n";
}

print "Type\tID\tChromosome\tStart\tEnd\tValue\n";	

while(<STDIN>)
{
	last if (/^ProbeSetName\tChromosome\tPosition\tLog2Ratio\tWeightedLog2Ratio\tSmoothSignal/);
}

my %probesets;
while(<STDIN>)
{
	last if  (/^#/);
	my ($ProbeSetName, $Chromosome, $Position, $Log2Ratio, $tWeightedLog2Ratio, $tSmoothSignal) = split(/\t/);
	$Chromosome = "X" if ($Chromosome eq "24");
	$Chromosome = "Y" if ($Chromosome eq "25");
	$probesets{$ProbeSetName} = "$Chromosome\t$Position";	
	next if ($Log2Ratio eq 'nan');
	next if (defined $chr and $Chromosome ne $chr);
	next if (defined $start and $Position < $start);
	next if (defined $end and $Position > $end);
	print "LRR\t$ProbeSetName\t$Chromosome\t$Position\t$Position\t$Log2Ratio\n";
}

while(<STDIN>)
{
	last if (/^SegmentID\tChromosome\tStartPosition\tStopPosition\tMarkerCount\tMeanMarkerDistance\tState\tConfidence/)
}

while(<STDIN>)
{
	last if (/^#/);
	my ($SegmentID, $Chromosome, $StartPosition, $StopPosition, $MarkerCount, $MeanMarkerDistance, $State, $Confidence) = split(/\t/);
	$Chromosome = "X" if ($Chromosome eq "24");
	$Chromosome = "Y" if ($Chromosome eq "25");
	next if (defined $chr and $Chromosome ne $chr);
	next if (defined $start and $StopPosition < $start);
	next if (defined $end and $StartPosition > $end);

	print "CN\t$SegmentID\t$Chromosome\t$StartPosition\t$StopPosition\t$State\n";	
}	

while(<STDIN>)
{
	last if (/^Index\tProbeSetName\tCall\tConfidence\tForcedCall\tASignal\tBSignal/)
}

while(<STDIN>)
{
	last if (/^#/);
	my ($Index, $ProbeSetName, $Call, $Confidence, $ForcedCall, $ASignal, $BSignal, $SignalStrength, $Contrast) = split(/\t/);
	die "ERROR: ProbeSetName not found: $_" if (!exists $probesets{$ProbeSetName});
	my ($Chromosome, $Position) = split(/\t/, $probesets{$ProbeSetName});
	next if (defined $chr and $Chromosome ne $chr);
	next if (defined $start and $Position < $start);
	next if (defined $end and $Position > $end);

	#my $baf = ($AllelePeaks0 - 250)/(750-250);
	my $baf = $ASignal / ($ASignal + $BSignal);
	$baf = 0 if ($baf < 0);
	$baf = 1 if ($baf > 1);
	print "BAF\t$ProbeSetName\t$Chromosome\t$Position\t$Position\t$baf\n";
}	

while(<STDIN>) {}; # read pipe to the end to avoid SIGPIPE error status 141
