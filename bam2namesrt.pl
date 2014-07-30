#!/usr/bin/perl -I /home/ksieber/scripts/ -I /home/ksieber/perl5/lib/perl5/
use warnings;
use strict;
use File::Basename;
use run_cmd;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
      'input=s',
      'output_prefix=s',
      'output_dir',
      'help',
);

if($options{help}){
   die "Help: This script will take a bam to sort based on names. 
      --input=           <BAM> to name sort.
        --sort_mem=		[500]
      --output_prefix=  prefix.bam
      --output_dir=     directory to put output
      --Qsub=			<0|1> [0] 1= Qsub the sort. 
        --sub_mem=		[6G] Needs to reflect changes in --sort_mem.
      --help\n";
}

my $input = $options{input} ? $options{input} : $ARGV[0];
my ($bam,$path,$suf)=fileparse($input,(".srt.bam",".sort.bam",".bam");
my $prefix = $options{output_prefix} ? $options{output_prefix} : $bam;
my $dir = $options{output_dir} ? $options{output_dir} : $path;
$dir =~ s/\/$//;
my $out = "$dir/$prefix";
my $sort_mem = $options{sort_mem} ? $options{sort_mem} : "foo";

my $cmd = "samtools sort -m $sort_mem";




