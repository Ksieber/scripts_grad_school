#!/usr/local/bin/perl -w

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
   	      'input=s',
	      'output_prefix=s',
	      'output_dir=s',
	      'h|help',
);

if ($options{help}){die "\nHELP: This script will use picard tools to calculat the insert size for a sorted.bam.
\t--input=            REQUIRED. Must be a sorted.bam
\t--output_prefix=    **Prefix**.insert.data
\t--output_dir=       Directory where to write the data if not the current working directory.";}

if (!$options{input}){die "ERROR: Must give a input.  Use: --input=some.sorted.bam\n";}
if (!$options{output_prefix}){
    $options{input}=~/(\w*)\.s.rt.*\.bam$/;
    $options{output_prefix}=$1;
}
my $output_prefix= $options{output_dir} ? "$options{output_dir}/$options{output_prefix}" : "$options{output_prefix}";
`java -Xmx3g -jar /home/jdhotopp/bin/Picard/picard-tools-1.48/CollectInsertSizeMetrics.jar AS=true I=$options{input} H=$output_prefix.std_insert.histogram O=$output_prefix.std_insert.metrics M=0 VALIDATION_STRINGENCY=SILENT DEVIATIONS=20`;
if ($? != 0 ){die "Fatal Error line 24\n";
`java -Xmx3g -jar /home/jdhotopp/bin/Picard/picard-tools-1.48/CollectInsertSizeMetrics.jar AS=true I=$options{input} H=$output_prefix.lrg_insert.histogram O=$output_prefix.lrg_insert.metrics M=0 VALIDATION_STRINGENCY=SILENT DEVIATIONS=1000000000000000000`;
if ($? != 0 ){die "Fatal Error line 26\n";
