#!/usr/local/bin/perl
use warnings;
use strict;
use run_cmd;
use setup_input;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
        'input=s',
        'input_list=s',
        'sort=s',
        'output_prefix=s',
        'output_dir=s',
        'Qsub=s',
        'region=s',
        'A=s',
        'd=s',
        'M_M=s',
        'M_UM=s',
        'help',
        );

if($options{help}){die
    "HELP: This script will take a bam file and calculate the coverage.
        --input=         	Input bam (Position sorted already)
        --input_list=     	List of bams to process.
        --sort=			<0|1> [0] 1= Sort by position. In order to calc. mpilupe input must be position sorted.
        --output_prefix=  	[filename] Prefix for the output file.
        --output_dir=    	Name the directory for output.  
        --region=         	chr1:100-200. Use to look @ reads only in this region. Highly recommended to use.
        --M_M=             	<0|1> [0] 1= calculate Mapped_Mapped reads only. 
        --M_UM=            	<0|1> [0] 1= calculate Mapped_UN-Mapped reads only. Works for LGT reads.
        --d=			[10000] Max Coverage per base. 
        --A=              	<0|1> [1] 1= Count anomalous read pairs (LGT).
        --Qsub=			<0|1> [0] 1= Qsub the mpileup foreach input.\n";
}
#### Check input
if(!$options{input} && !$options{input_list}){die "ERROR: Must give an input bam with --input=<BAM>. Try again.\n";}
if($options{M_M} && $options{M_UM}){die "ERROR: Can only use one: M_M or M_UM at once. Try again.\n";}

## Setup options
my $qsub = $options{Qsub} ? $options{Qsub} : 0;
my $sort = $options{sort} ? $options{sort} : 0;
my $MM = $options{MM} ? $options{MM} : 0;
my $MU = $options{MU} ? $options{MU} : 0;
my $A = $options{A} ? $options{A} : 1;
my $d = $options{d} ? $options{d} : 100000;
my $view = "-hu";              ## Default
my $mpileup = "-Ad $d";    	   ## Default
if($MM==1){
    $view = "-huF0xC";
    $mpileup = "-d $d";
}
if($MU==1){
    $view = "-huf0x8 -F0x4";
}  
if($A==0){
    $mpileup = "-d $d";
}


my $input = setup_input(\%options);
run_cmd("mkdir -p $options{output_dir}");

foreach my $bam (@$input){
	my ($filename,$directory) = fileparse($bam, ('.srt.bam','.bam'));
    my $prefix = $options{output_prefix} ? $options{output_prefix} : $filename;
    my $dir = $options{output_dir} ? $options{output_dir} : $directory;
    my $out = "$dir/$prefix";
    my $cmd;
    if($sort==1){
    	$cmd = "samtools sort -m 5000000000 -o $bam - | samtools view $view - \'$options{region}\' | samtools mpileup $mpileup - > $out\.mpileup"
    } else {
    	$cmd = "samtools view $view $bam \'$options{region}\' | samtools mpileup $mpileup - > $out\.mpileup";
 	}
    if($qsub==1){Qsub($cmd); next;} 
    else {run_cmd($cmd);}
}	
