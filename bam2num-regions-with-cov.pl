#!/usr/bin/perl -I /home/ksieber/scripts/ -I /home/ksieber/perl5/lib/perl5/
use warnings;
use strict;
use run_cmd;
use File::Basename;
use setup_input;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'input_list|I=s', 'min_cov|c=i', 'window_size|w=i', 'step_size|s=i', 'sort|S=i', 'threads|t=i', 'sort_mem=s', 'ref|r=s', 'help|?' )
    or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} ) {
    die "This script will calculate the number of regions with an average >= --min_cov over --window_size.
	--input|i=          
        --input_list|I=     
        --min_cov|c=		[2]
	--window_size|w= 	[25]
	--step_size|s= 		[10]
	--sort|S=		[1] 
	--ref|r=	
	--threads|t=		[1]
	--sort_mem=		[1G]
	--help|?
	\n";
}
if ( !$options{input} && !$options{input_list} ) { die "Error: Must pass an input with --input or --input_list. Please try again.\n"; }

my $min_cov  = defined $options{min_cov}     ? $options{min_cov}     : 2;
my $window   = defined $options{window_size} ? $options{window_size} : 25;
my $step     = defined $options{step_size}   ? $options{step_size}   : 10;
my $sort     = defined $options{sort}        ? $options{sort}        : 1;
my $threads  = defined $options{threads}     ? $options{threads}     : 1;
my $sort_mem = defined $options{sort_mem}    ? $options{sort_mem}    : "1G";

my $inputs = setup_input( \%options );
foreach my $input (@$inputs) {
    my ( $fn, $dir, $suf ) = fileparse( $input, ( '.psort.bam', '.srt.bam', '.bam' ) );
    if ( $sort == 1 ) {
        my $cmd1 = "/home/ksieber/scripts/bam2srtbai.pl --input=$input --threads=$threads --sort_mem=$sort_mem";
        run_cmd($cmd1);
    }
    my $cmd2 = ( $sort == 1 ) ? "/home/ksieber/scripts/bam2mpileup.pl $dir/$fn\.psort.bam" : "/home/ksieber/scripts/bam2mpileup.pl $input";
    run_cmd($cmd2);
    my @regions
        = ( $sort == 1 )
        ? `/home/ksieber/scripts/sliding_window_coverage_v2.pl --input=$dir/$fn\.psort.mpileup --window_size=$window --min_cov=$min_cov --step_size=$step`
        : `/home/ksieber/scripts/sliding_window_coverage_v2.pl --input=$dir/$fn\.mpileup --window_size=$window --min_cov=$min_cov --step_size=$step`;
    chomp(@regions);
    my $number_of_windows = scalar(@regions);
    print STDOUT "$fn$suf:\t$number_of_windows\tregions of: $window bp with: $min_cov coverage.\n";
}
