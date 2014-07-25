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

# read kgXref, knownCanonical to determine UCSC canonical transcripts affected by variant
my %kgID2refSeq;
open(G,"$ENV{HOME}/generic/data/hg19/hg19.kgXref.txt") or die "could not open file $ENV{HOME}/generic/data/hg19/hg19.kgXref.txt";
while(<G>)
{
	chomp;
	my ($kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refSeq, $protAcc, $description, $rfamAcc, $tRnaName) = split(/\t/);

	$kgID2refSeq{$kgID} = $refSeq if ($refSeq);
}
close(G);
INFO(scalar(keys(%kgID2refSeq))." gene descriptions read from file $ENV{HOME}/generic/data/hg19/hg19.kgXref.txt");

my %canonical;
open(G,"$ENV{HOME}/generic/data/hg19/hg19.knownCanonical.txt") or die "could not open file $ENV{HOME}/generic/data/hg19/hg19.knownCanonical.txt";
<G>; # skip header
while(<G>)
{
	chomp;
	my ($chrom, $chromStart, $chromEnd, $clusterId, $transcript, $protein) = split(/\t/);
	
	$canonical{$kgID2refSeq{$transcript}} = 1 if ($kgID2refSeq{$transcript});
}
close(G);
INFO(scalar(keys(%canonical))." canonical genes read from file $ENV{HOME}/generic/data/hg19/hg19.knownCanonical.txt");

INFO("Processing file $vcf_in...");

my $vcf = Vcf->new(file => "$vcf_in");
$vcf->parse_header();

my (@samples) = $vcf->get_samples();

print "sample\tchr\tpos\tdbSNP\tref\talt\tgene\timpact\teffect\tdp_rem\tdp_leu\taf\n";
while (my $line = $vcf->next_line())
{
	my $x = $vcf->next_data_hash($line);
	my ($gene, $add_genes, $impact, $effect, $affected_exon, $aa_change) = get_impact($x->{INFO}{EFF});

	foreach my $s (@samples) 
	{
		my ($ad_ref, $ad_alt) = split(",", $x->{gtypes}{$s}{AD});
		$ad_alt = $ad_ref if (!defined $ad_alt);

		next if ($rem_samples{$s} || $ad_alt < 6);

		print "$s\t".$x->{CHROM}."\t".$x->{POS}."\t".$x->{ID}."\t".$x->{REF}."\t".$x->{ALT}->[0]."\t$gene\t$impact\t$effect\t$ad_ref\t$ad_alt\t".($ad_alt/($ad_ref+$ad_alt))."\n"	
	}
}

# ------------------------------------------

sub get_impact
{
	my $effs = shift or die "ERROR: effect not specified";

	# determine all genes impacted by variants
	my (%genes_by_impact, %all_genes, $combined_effect, $combined_impact, %affected_exons, %aa_changes);
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
			$affected_exons{$gene_name}{$exon}{$transcript} = 1;
			if ($canonical{$transcript})
			{
				$affected_exons{$gene_name}{'canonical'}{$exon}{$transcript} = 1;
			}
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

	my @aff_exons;
	foreach my $g (keys(%affected_exons))
	{
		if (exists $affected_exons{$g}{'canonical'}) # known canonical transcript for this gene?
		{
			foreach my $e (keys(%{$affected_exons{$g}{'canonical'}}))
			{
				my @transcripts;
				foreach my $t (keys(%{$affected_exons{$g}{'canonical'}{$e}}))
				{
					push(@transcripts, "$g:$t");
				}
				push(@aff_exons, "$e (".join(";", @transcripts).")");
			}
		}
		else
		{
			foreach my $e (keys(%{$affected_exons{$g}}))
			{
				next if ($e eq 'canonical');

				my @transcripts;
				foreach my $t (keys(%{$affected_exons{$g}{$e}}))
				{
					push(@transcripts, "$g:$t");
				}
				push(@aff_exons, "$e (".join(";", @transcripts).")");
			}
			
		}
	}

	return ($gene, $add_genes, $combined_impact, $combined_effect, 
			@aff_exons > 0 ? join(",", @aff_exons) : "", join(";", keys(%aa_changes)));
}
	