use warnings FATAL => qw( all );
use strict;

use Getopt::Long;
use File::Basename;

my ($template, $data_dir, $file_pattern, $gfx_output_file);
GetOptions
(
	"template=s" => \$template,  # template configuration file
	"data-dir=s" => \$data_dir,  # directory containing circos data files
	"file-pattern=s" => \$file_pattern,  # regular expression matching data file names
	"gfx-output-file=s" => \$gfx_output_file  # graphics output file
);

die "ERROR: --template not specified\n" if (!$template);
die "ERROR: --data-dir not specified\n" if (!$data_dir);
die "ERROR: --file-pattern not specified\n" if (!$file_pattern);
die "ERROR: --gfx-output-file not specified\n" if (!$gfx_output_file);

opendir(DIR, $data_dir) or die "ERROR: Could not read directory $data_dir\n";

my @files;
while(my $file = readdir(DIR))
{
	next unless ($file =~ /$file_pattern/);
	push(@files, "$data_dir/$file");
}
closedir(DIR);

print STDERR "Found ".scalar(@files)." files in directory $data_dir\n";

print STDERR "Processing template $template...\n";

my $gfx_output_file_dir = dirname($gfx_output_file);
my $gfx_output_file_name = basename($gfx_output_file);

open(T, "$template") or die "ERROR: Could not open template configuration file $template\n";
while(<T>)
{
	if (/\s*#/)
	{
		print $_;
		next;
	}
	
	s/\[\[OUT_DIR\]\]/$gfx_output_file_dir/;
	s/\[\[OUT_FILE\]\]/$gfx_output_file_name/;
	
	if (/\[\[PLOTS\]\]/)
	{
		my $r = 0.97;
		foreach my $f (@files)
		{
			print "\t<plot>\n";
			print "\t\ttype = highlight\n";
			print "\t\tcolor = white\n";
			print "\t\tfile = $f\n";
        	print "\t\tr0 = ",$r-0.026,"r\n";
        	print "\t\tr1 = ",$r,"r\n";
			print "\t</plot>\n";
			
			$r = $r - 0.03;
		}
		next;
	}

	if (/\[\[HIGHLIGHTS\]\]/)
	{
		my $r = 0.99;
		foreach my $f (@files)
		{
			print "\t<highlight>\n";
			print "\t\tfile = $f\n";
        	print "\t\tr0 = ",$r-0.02,"r\n";
        	print "\t\tr1 = ",$r,"r\n";
			print "\t</highlight>\n";
			
			$r = $r - 0.02;
		}
		next;
	}
	
	print $_;
}
close(T);