use warnings FATAL => qw( all );
use strict;

use Getopt::Long;
use File::Basename;

my ($template, $data_dir, $order_file, $file_pattern, $gfx_output_file);
GetOptions
(
	"template=s" => \$template,  # template configuration file
	"data-dir=s" => \$data_dir,  # directory containing circos data files
	"order-file=s" => \$order_file, # file with order of samples (outward to inward)
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
	push(@files, "$file");
}
closedir(DIR);

# order samples
if ($order_file)
{
	my @ordered_files;
	open(O, "$order_file") or die "ERROR: Could not read sample order file $order_file\n";
	while(<O>) 
	{
		chomp; next if (!$_ or /^\s*#/);
		
		if ($_ eq "---") 
		{
			print STDERR "Adding empty track (separator) to list of ordered samples.\n";
			push(@ordered_files, $_);
			next;
		}
		
		my $found = 0;
		for (my $i = 0; $i < @files; $i ++)
		{
			if ($files[$i] =~ /$_/) 
			{
				print STDERR "Adding $files[$i] to list of ordered samples.\n";
				push(@ordered_files, $files[$i]);
				splice(@files, $i, 1);
				$found = 1;
				last;
			}
		}
		
		die("ERROR: Sample identifier $_ not found among Circos files present in directory $data_dir\n")
			if (!$found); 
	}
	close(O);
	
	push(@ordered_files, @files);
	@files = @ordered_files;
}

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
			if ($f ne "---")
			{
				print "\t<plot>\n";
				print "\t\ttype = highlight\n";
				print "\t\tcolor = white\n";
				print "\t\tfile = $data_dir/$f\n";
	        	print "\t\tr0 = ",$r-0.011,"r\n";
	        	print "\t\tr1 = ",$r,"r\n";
				print "\t</plot>\n";				
			}
			
			$r = $r - 0.013;
		}
		next;
	}

	if (/\[\[HIGHLIGHTS\]\]/)
	{
		my $r = 0.99;
		foreach my $f (@files)
		{
			if ($f ne "---")
			{
				print "\t<highlight>\n";
				print "\t\tfile = $data_dir/$f\n";
	        	print "\t\tr0 = ",$r-0.009,"r\n";
	        	print "\t\tr1 = ",$r,"r\n";
				print "\t</highlight>\n";
			}
			
			$r = $r - 0.01;
		}
		next;
	}
	
	print $_;
}
close(T);