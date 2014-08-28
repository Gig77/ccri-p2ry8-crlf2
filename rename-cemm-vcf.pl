use warnings FATAL => qw( all );
use strict;

my $dir = "/home/STANNANET/christian.frech/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_zahyee0Zooxiphoogo5ooMee3mai8uy8/";
opendir(DIR, "$dir") or die $!;

while (my $file = readdir(DIR)) {

	# We only want files
	next unless (-f "$dir/$file");

	# Use a regular expression to find files
	if ($file =~ /variant_calling_somatic_(.+)__(.+)_annotated.vcf$/)
	{
		my ($s1, $s2) = ($1, $2);
		my $target;
		
		if ($s1 =~ /^.+(C|_Remission)$/ and $s2 =~ /^(.+)(D|_Diagnosis)$/)
		{
			$target = "/home/STANNANET/christian.frech/p2ry8-crlf2/data/mutect/$1_rem_dia.somatic.vcf";
		}
		elsif ($s1 =~ /^.+(C|_Remission)$/ and $s2 =~ /^(.+)(R|R1|_Relapse)$/)
		{
			$target = "/home/STANNANET/christian.frech/p2ry8-crlf2/data/mutect/$1_rem_rel.somatic.vcf";
		}
		elsif ($s1 =~ /^.+(C|_Remission)$/ and $s2 =~ /^(.+)R2$/)
		{
			$target = "/home/STANNANET/christian.frech/p2ry8-crlf2/data/mutect/$1_rem_rel2.somatic.vcf";
		}
		else
		{
			die "Could not parse filename: $file\n";
		}
		
		system("rm -f $target");
		my $cmd = "ln -s $dir/$file $target";
		print("$cmd\n");
		system($cmd);
	}
}
closedir(DIR);