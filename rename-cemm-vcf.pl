use warnings FATAL => qw( all );
use strict;

system("ln -s /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0024_HD_BCP_ALL_xiewahYoojiev7ahyoh4ohjiseiKepoh/variant_calling_somatic_545_Remission__m1963_545_Relapse_annotated.vcf /mnt/projects/p2ry8-crlf2/data/mutect/m1963-545-rel_rem_xeno.somatic.vcf");
system("ln -s /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0024_HD_BCP_ALL_xiewahYoojiev7ahyoh4ohjiseiKepoh/variant_calling_somatic_G_Remission__m1977_G_Dx_Diagnosis_annotated.vcf /mnt/projects/p2ry8-crlf2/data/mutect/m1977-G-dia_rem_xeno.somatic.vcf");
system("ln -s /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0024_HD_BCP_ALL_xiewahYoojiev7ahyoh4ohjiseiKepoh/variant_calling_somatic_Y3767_Remission__m1967_Y_Relapse_annotated.vcf /mnt/projects/p2ry8-crlf2/data/mutect/m1967-Y-rel_rem_xeno.somatic.vcf");
system("ln -s /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0024_HD_BCP_ALL_xiewahYoojiev7ahyoh4ohjiseiKepoh/variant_calling_somatic_715_Remission__M1957_715_Relapse1_annotated.vcf /mnt/projects/p2ry8-crlf2/data/mutect/m1957-715-rel_rem_xeno.somatic.vcf");
system("ln -s /mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0024_HD_BCP_ALL_xiewahYoojiev7ahyoh4ohjiseiKepoh/variant_calling_somatic_545_Remission__m1964_545Rn_Relapse_annotated.vcf /mnt/projects/p2ry8-crlf2/data/mutect/m1964-545-rel_rem_xeno.somatic.vcf");

my $dir = "/mnt/projects/p2ry8-crlf2/data/mutect/www.biomedical-sequencing.at/projects/BSA_0026_P2RY8_CRLF2_ALL_zahyee0Zooxiphoogo5ooMee3mai8uy8";
opendir(DIR, "$dir") or die $!;

while (my $file = readdir(DIR)) {

	# We only want files
	next unless (-f "$dir/$file");

	# Use a regular expression to find files
	if ($file =~ /variant_calling_somatic_(.+)__(.+)_(snpeff|annotated).vcf$/)
	{
		my ($s1, $s2) = ($1, $2);
		print "$s1 $s2\n";
		my $target;

		if ($s1 eq "842_Remission" and $s2 eq "841_Diagnosis")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/841_rem_dia.somatic.vcf";			
		}
		elsif ($s1 eq "DL2_Remission" and $s2 eq "19981_DL2_R_Relapse")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/DL2_rem_rel.somatic.vcf";			
		}
		elsif ($s1 eq "MA5_Remission" and $s2 eq "19319_MA5_R_Relapse")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/MA5_rem_rel.somatic.vcf";			
		}
		elsif ($s1 eq "GI8_Remission" and $s2 eq "19551_GI8_R_Relapse")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/GI8_rem_rel.somatic.vcf";			
		}
		elsif ($s1 eq "BJ17183_Remission" and $s2 eq "14367_Relapse")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/BJ17183_rem_rel.somatic.vcf";			
		}
		elsif ($s1 eq "KE17247_Remission" and $s2 eq "15721_Relapse")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/KE17247_rem_rel.somatic.vcf";			
		}
		elsif ($s1 eq "DS10898_Remission" and $s2 eq "5143_Relapse")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/DS10898_rem_rel.somatic.vcf";			
		}
		elsif ($s1 eq "VS14645_Remission" and $s2 eq "9931_Relapse")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/VS14645_rem_rel.somatic.vcf";			
		}
		elsif ($s1 eq "SE1528_5_Remission" and $s2 eq "13977_Relapse")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/SE15285_rem_rel.somatic.vcf";			
		}
		elsif ($s1 eq "AL9890_Remission" and $s2 eq "AL9890_Diagnosis")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/AL9890_rem_dia.somatic.vcf";			
		}
		elsif ($s1 eq "AL9890_Remission" and $s2 eq "AL9890_Relapse")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/AL9890_rem_rel.somatic.vcf";			
		}
		elsif ($s1 eq "GL11356_Remisson" and $s2 eq "GL11356_Diagnosis")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/GL11356_rem_dia.somatic.vcf";			
		}
		elsif ($s1 eq "GL11356_Remisson" and $s2 eq "GL11356_Relapse")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/GL11356_rem_rel.somatic.vcf";			
		}
		elsif ($s1 =~ /^.+(C|_Remission)$/ and $s2 =~ /^(.+)(D|_Diagnosis)$/)
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/$1_rem_dia.somatic.vcf";
		}
		elsif ($s1 =~ /^.+(C|_Remission)$/ and $s2 =~ /^(.+)(R|R1|_Relapse)$/)
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/$1_rem_rel.somatic.vcf";
		}
		elsif ($s1 =~ /^.+(C|_Remission)$/ and $s2 =~ /^(.+)R2$/)
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/$1_rem_rel2.somatic.vcf";
		}
		elsif ($s1 =~ /^.+(C|_Remission)$/ and $s2 =~ /^(.+)R3$/)
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/$1_rem_rel3.somatic.vcf";
		}
		elsif ($s1 eq "737C" and $s2 eq "737_Relapse3")
		{
			$target = "/mnt/projects/p2ry8-crlf2/data/mutect/737_rem_rel3.somatic.vcf";
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