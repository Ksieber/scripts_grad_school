#!/usr/local/bin/perl

## This is not complete! Needs work
## I want to make it smart to be able to handle multiple chrs, and print windows when large gaps in cov. 

use warnings;
no warnings 'uninitialized';
use strict;
use run_cmd;
use Carp;
$Carp::MaxArgLen = 0;
use File::Basename;
use lib qw(/local/projects-t3/HLGT/scripts/lgtseek/lib/);
use LGTSeek;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
        'input|i=s',
        'sort|s=i',
        'window_size|w=i',
        'step_size|p=i',
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

if($options{help}){help();}
if(!$options{input}){die "ERROR: Must give an input file with: --input= BAM or mpileup file."}
unless(-e $options{input}){die "ERROR: Input file can't be found.";}
if(!$options{window_size} || !$options{step_size}){die "ERROR: Must pass the --window_size & --step_size. Use --Help for more info.\n";}

my $lgtseq = LGTSeek->new2(\%options);
our $window_size = defined $options{window_size} ? $options{window_size} : "10";
our $step_size = defined $options{step_size} ? $options{step_size}-1 : "5";

my $sort = defined $options{sort} ? "$options{sort}" : "0";
my $qsub = defined $options{Qsub} ? "$options{Qsub}" : "0";
my $region = defined $options{region} ? "'$options{region}'" : undef;
my $A = defined $options{A} ? $options{A} : "1";
my $d = defined $options{d} ? $options{d} : 100000;
my $view = "-hu";              ## Default
my $mpileup = "-Ad $d";    	   ## Default
if($A==0){ $mpileup = "-d $d"; }

my $window;					# keys: coverage (total coverage seen in window), start (actual start), stop (expected stop), bases (total in window seen), chr
my $next_window;
my $fh = open_input($options{input});

while(<$fh>){
	chomp;
	my $line = $_; 
	my ($chr,$pos,$coverage)=(split/\t/,$line)[0,1,3];
	if($chr ne $window->{'chr'}){
		if(defined $window->{'chr'}){print calc_cov($window)."\n";}
		$window = init_window($line);
		next;
	}
	if($pos==($window->{'start'}+$step_size)){$next_window = init_window($line);}
	if($pos > $next_window->{'start'} && $pos < $next_window->{'stop'}){
		$next_window->{'coverage'}+=$coverage;
		$next_window->{'bases'}++;
	}
	# if within window size, add to window
	if($pos > $window->{'start'} && $pos < $window->{'stop'}){
		$window->{'coverage'}+=$coverage;
		$window->{'bases'}++;
		next;
	}
	# if outside of window size, print window and move to next_window.
	if($pos >= $window->{'stop'}){
		print calc_cov($window)."\n";
		if($pos>$next_window->{'stop'}){$next_window = undef;}
		$window = $next_window;
		$next_window = init_window($line);
	}
}

print calc_cov($window)."\n";
close $fh;





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
	my $window = shift;
	my $average = sprintf("%.2f", $window->{'coverage'} / $window->{'bases'});
	my $actual_stop = ($window->{'start'}+($window->{'bases'}-1));
	my $step_pos = $window->{'start'} + ($step_size-1);
	return "$average\t$window->{chr}\:$window->{start}\-$actual_stop\t$window->{chr}\t$step_pos";
}

sub init_window {
	my $line = shift;
	my ($chr,$pos,$coverage)=(split(/\s/,$line))[0,1,3];
	my %window;
	$window{chr}=$chr;
	$window{start}=$pos;
	$window{stop}=$pos+$window_size;
	$window{coverage}=$coverage;
	$window{bases}=1;
	return \%window;
}

sub help {
	die 
	"Help: This script will calculate coverage with a sliding window. 
		--input=	Mpileup or bam.
		--help\n";
}


