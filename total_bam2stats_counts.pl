#!/usr/bin/perl -I /home/ksieber/scripts/ -I /home/ksieber/perl5/lib/perl5/
use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'help|?' ) or die "Error: Unrecognized command line option. Please try again.\n";

my @columns = ( "bam", "total", "MM", "MU", "UU", "SC" );
my %counts;
map { $counts{$_} = 0; } (@columns);
if ( $options{help} ) { &help; }
if ( !$options{input} && !$ARGV[0] ) { die "Error: Must pass an input bam2stats.txt file with --input or $ARGV[0]\. Please try again.\n"; }
my $input = defined $options{input} ? $options{input} : $ARGV[0];

open( IN, "<", "$input" ) or die "Error: Unable to open input: $options{input}\n";
while (<IN>) {
    chomp;
    next if($_=~/^BAM/);
    next if($_=~/^Total-Bams/);
    my ( $bam, $total, $MM, $MU, $UU, $SC ) = split;
    $counts{bam}++;
    $counts{'total'} = $counts{'total'} + $total;
    $counts{'MM'}    = $counts{'MM'} + $MM;
    $counts{'MU'}    = $counts{'MU'} + $MU;
    $counts{'UU'}    = $counts{'UU'} + $UU;
    $counts{'SC'}    = $counts{'SC'} + $SC;
}
close IN;

my $OUT = *STDOUT;
print $OUT "Total-Bams:";
map { print $OUT "$counts{$_}\t"; } (@columns);
print $OUT "\n";

sub help {
    die "This script will calculate the total numbers of a bam2stats.txt file. Output => STDOUT.
	--input|i= 	bam2stats.txt file.
	--help|?\n";
}
