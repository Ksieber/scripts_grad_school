#!/usr/local/bin/perl

## This is not complete! Needs work
## I want to make it smart to be able to handle multiple chrs, and print windows when large gaps in cov. 

use warnings;
use strict;
use run_cmd;
use Carp;
$Carp::MaxArgLen = 0;
use File::Basename;
use lib qw(/local/projects-t3/HLGT/scripts/lgtseek/lib/);
use LGTSeek;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
        'input=s',
        'sort=i',
        'window_size=i',
        'step_size=i',
        'output_prefix=s',
        'output_dir=s',
        'Qsub=i',
        'region=s',
        'A=i',
        'd=i',
        'M_M=i',
        'M_UM=i',
        'help',
        );

if($options{help}){&help;}
if(!$options{input}){confess "ERROR: Must give an input file with: --input= BAM or mpileup file."}
unless(-e $options{input}){confess "ERROR: Input file can't be found.";}

my $lgtseq = LGTSeek->new2(\%options);
my $window_size = $options{window_size};
my $step_size = $options{step_size}-1;

my $sort = defined $options{sort} ? "$options{sort}" : "0";
my $qsub = defined $options{Qsub} ? "$options{Qsub}" : "0";
my $region = defined $options{region} ? "'$options{region}'" : undef;
my $A = defined $options{A} ? $options{A} : "1";
my $d = defined $options{d} ? $options{d} : 100000;
my $view = "-hu";              ## Default
my $mpileup = "-Ad $d";    	   ## Default
if($A==0){ $mpileup = "-d $d"; }

my $fh = &open_input($options{input});

my %window;				## $window{chr}{$position}=$coverage;
my $previous_chr;
my $window_start_position = 0;
my @window_positions;


while(<$fh>){
	chomp; 
	my ($chr,$pos,$coverage)=(split)[0,1,3];
	
	if($chr ne '$previous_chr'){
		print "&calc_cov(\%window)";
		$window{$previous_chr}=undef;
	}

	if($pos < ($window_start_position + $window_size)){
		$window{$chr}{$pos}=$coverage;
		push(@window_positions,$pos);
	}

	if($pos == ($window_start_position + $window_size)){
		print "&calc_cov(\%window)\n";
		$window_start_position = $window_start_position + $step_size;

	}

	if($pos > ($window_start_position + $window_size)){
		print "&calc_cov(\%window)\n";
		$window_start_position = $pos; 

	}
}



sub open_input{
	my $raw_input = shift;
	my ($fn,$path,$suffix)=fileparse($raw_input,(@{$lgtseq->{mpileup_suffix_list}},@{$lgtseq->{bam_suffix_list}}));
	
	if($lgtseq->empty_chk({input => $raw_input})==1){confess "ERROR: Input file:$raw_input is empty.";}
	## Filehandle to return.
	my $fh;
	## Mpileup input
	map {
		if($suffix eq $_){
			open($fh, "<","$raw_input") || confess "ERROR: Can't open input: $raw_input because: $!\n"; 
			return $fh;
		}
	} @{$lgtseq->{mpileup_suffix_list}};

	## Position sorted bam
	my @psort_bams = ('_pos-sort.bam','_psort.bam','srt.bam');
	map {
		if($suffix eq $_){
			open($fh, "-|","samtools view $view $raw_input $region | samtools mpileup $mpileup - ") || confess "ERROR: Can't open input: $raw_input because: $!\n"; 
			return $fh;
		}
	} @psort_bams;


	## Name sorted bam || --sort=1
	map {
		if($suffix eq $_ || $sort==1){
			print STDERR "Sorting input bam . . .\n";
			open($fh, "-|","samtools sort -m 5G -o $raw_input - | samtools view $view - $region | samtools mpileup $mpileup - ") || confess "ERROR: Can't open input: $raw_input because: $!\n"; 
			return $fh;
		}
	}@{$lgtseq->{bam_suffix_list}};
}

sub calc_cov {
	my %window = shift;
	foreach my $chr (keys %window){
		foreach $pos (keys %{$window{$chr}}){
			
		}
	}
}

sub help {
	die 
	"Help: This script will calculate coverage with a sliding window. 
		--input=	Mpileup or bam.
		--help\n";
}

