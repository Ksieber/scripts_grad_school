#!/usr/bin/perl
use warnings;
use strict;
use lib "/local/projects-t3/HLGT/scripts/lgtseek/lib/";      ### May need to change this depending on where the script is being run
use LGTSeek;
use Bio::DB::Fasta;
use Bio::Graphics;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use File::Basename;
use run_cmd;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
	'ref1=s',
	'ref2=s',
	'region1=s',
	'region2=s',
	'bam1=s',
	'bam2=s',
	'output_dir=s',
	'output_prefix=s',
	'help'
);

if($options{help}){
	die 
"Help: This script will take a 2 bams mapped @ different references and merge the reference & map at the merged reference. This is useful for merging LGT regions.
		--ref1=				Reference 1 fasta.
		--ref2=				Reference 2 fasta.
		--region1=			<chr#:100-200>
		--region2=			<chr#:100-200>
		--bam1=				bam1. 5' side of merger. Assumes position sorted.
		--bam2=				bam2. 3' side of merger. Assumes position sorted.
		--sort1=			<0|1> [0] 1= Position Sort bam1.
		--sort2=			<0|1> [0] 1= Position Sort bam2.
		--output_dir=		Directory for output. 		[bam1 dir]
		--output_prefix=	Prefix for the output fasta & bam. 		[bam1-merged-bam2]
";}

if(!$bam1 || !$bam2	|| !$ref1 || !$ref2){
	die "ERROR: Must pass 2 input references and bams. --bam1= --bam2= --ref1= --ref2=\n";
}

my ($fn1,$path1,$suff1) = fileparse($options{bam1},(".srt.bam",".bam"));
my ($fn2,$path2,$suff2) = fileparse($options{bam2},(".srt.bam",".bam"));

my $sort1 = defined $options{sort1} : "$options{sort1}" : "0";
my $sort2 = defined $options{sort2} : "$options{sort2}" : "0";

if($sort1==1){
	run_cmd("samtools sort $options{bam1} $path1/$fn1\_psort");
	run_cmd("samtools index $path1/$fn1\_psort.bam $path1/$fn1\_psort.bai")
	$options{bam1} = "$path1/$fn1\_psort.bam";
}
if($sort2==1){
	run_cmd("samtools sort $options{bam2} $path2/$fn2\_psort");
	run_cmd("samtools index $path2/$fn2\_psort.bam $path2/$fn1\_psort.bai")
	$options{bam1} = "$path1/$fn2\_psort.bam";
}

my $ref1 = Bio::DB::Fasta->new("$options{ref1}");
my $ref2 = Bio::DB::Fasta->new("$options{ref2}");

my $samtools-view1 = defined $options{region1} : "samtools view $options{bam1} $options{region1}" : "samtools view $options{bam1}";
my $samtools-view2 = defined $options{region2} : "samtools view $options{bam2} $options{region2}" : "samtools view $options{bam2}";
open(my $in1, "$samtools-view1 |") || die "ERROR: Can't open bam1: $options{bam1} because: $!\n";
open(my $in2, "$samtools-view2 |") || die "ERROR: Can't open bam1: $options{bam2} because: $!\n";

while(<$in1>){
	chomp; 
}

