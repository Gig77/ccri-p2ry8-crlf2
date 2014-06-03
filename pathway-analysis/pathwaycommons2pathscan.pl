use warnings FATAL => qw( all );
use strict;

use Log::Log4perl qw(:easy);

# read id mapping
my %id2sym;
open(M, "$ENV{HOME}/hdall/results/id-mappings.tsv") or die "ERROR: could not read id mappings\n";
while(<M>)
{
	chomp;
	my ($sym, $id) = split(/\t/);
	$id2sym{$id} = $sym;
}
close(M);
INFO(scalar(keys(%id2sym))." id mappgins read from file $ENV{HOME}/hdall/results/id-mappings.tsv");

# read kgXref
my (%uni2sym);
open(G,"$ENV{HOME}/generic/data/hg19/hg19.kgXref.txt") or die "could not open file $ENV{HOME}/generic/data/hg19/hg19.kgXref.txt";
while(<G>)
{
	chomp;
	my ($kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refSeq, $protAcc, $description, $rfamAcc, $tRnaName) = split(/\t/);
	next if (!$spID);

	$geneSymbol = $id2sym{$geneSymbol} if ($id2sym{$geneSymbol});
	$geneSymbol = $id2sym{$geneSymbol} if ($id2sym{$geneSymbol});
	
	$uni2sym{$spID} = $geneSymbol if ($spID);
}
close(G);
INFO(scalar(keys(%uni2sym))." id mappgins read from file $ENV{HOME}/generic/data/hg19/hg19.kgXref.txt");
# read biomart id mapping to get entrez ids and additional uniprot ids
my %sym2entrez;
open(G, "$ENV{HOME}/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv") or die "ERROR: could not read gene list\n";
while(<G>)
{
	chomp;
	my ($approved_symbol, $entrez_gene_id, $accession_numbers, $approved_name, 
		$previous_symbols, $previous_names, $aliases, $name_aliases, $entrez_gene_id_ncbi, $ensembl_id, $uniprot, $refseq, $ucsc) = split("\t");

	$approved_symbol = $id2sym{$approved_symbol} if ($id2sym{$approved_symbol});
	$uni2sym{$uniprot} = $approved_symbol if ($uniprot);

	my @symbols;
	push(@symbols, $approved_symbol) if ($approved_symbol);
	push(@symbols, split(", ", $previous_symbols)) if ($previous_symbols);
	push(@symbols, split(", ", $aliases)) if ($aliases);
	
	foreach my $s (@symbols)
	{
		$s = $id2sym{$s} if ($id2sym{$s});
		$sym2entrez{$s} = $entrez_gene_id if ($entrez_gene_id);
		$sym2entrez{$s} = $entrez_gene_id_ncbi if (!$entrez_gene_id and $entrez_gene_id_ncbi);
	}
}
close(G);

my ($read, $mapped, $lines) = (0, 0, 0);
my %pathways;
while(<>)
{
	chomp;

	#10090:  5HT4 type receptor mediated signaling pathway	datasource: panther; organism: 10090; id type: uniprot	P51829	P97288	Q3V1Q3
	my ($id, $description, $datasource, $organism, $id_type, $uniprot_ids) = /^([^:]+):\s+(.*?)\s+datasource: ([^;]+); organism: ([^;]+); id type: (\S+)\s+(.*)/;
	
	$lines ++;
	next if (!$id or $organism ne "9606");

	$id = "$datasource|$description";
	foreach my $u (split(/\s+/, $uniprot_ids))
	{
		$read ++;
		next if (!$uni2sym{$u}); 
	
		$pathways{$id}{$u} = $uni2sym{$u};
		$mapped ++;
	}
}
print STDERR "$lines input lines mapped to ".scalar(keys(%pathways)). " pathways\n";
print STDERR "Mapped $mapped out of $read UniProt IDs\n";

my $id = 1;
foreach my $p (keys(%pathways))
{
	my ($datasource, $description) = split('\|', $p);
	
	# MuSiC required format:
	# This is a tab-delimited file prepared from a pathway database (such as KEGG), with the columns: 
	# [path_id, path_name, class, gene_line, diseases, drugs, description] The latter three columns 
	# are optional (but are available on KEGG). The gene_line contains the "entrez_id:gene_name" of 
	# all genes involved in this pathway, each separated by a "|" symbol.
	# Example:
	# hsa00061      Fatty acid biosynthesis	      Lipid Metabolism     31:ACACA|32:ACACB|27349:MCAT|2194:FASN|54995:OXSM|55301:OLAH	
	print "PC4_$id\t";
	print "$datasource\t";
	print "$description\t";

	my @genes = values(%{$pathways{$p}});
	for (my $i = 0; $i < @genes; $i++)
	{
		my $entrezid = $sym2entrez{$genes[$i]};
		$entrezid = "0000" if (!$entrezid);
		
		print "$entrezid:$genes[$i]";
		print "|" if ($i < @genes-1);
	}
		
	print "\t\t\t\n";
	
	$id ++;
}