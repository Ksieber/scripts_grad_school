#!/usr/bin/perl
use strict;
use warnings;
no warnings 'uninitialized';

if ( !@ARGV ) { &help; }

use Getopt::Long qw (:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'lca_file=s', 'bam_file=s', 'good_list=s', 'bad_list=s', 'help|?' )
    or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} ) { &help; }

if ( !$options{lca_file} ) {
    die "ERROR: You must give an LCA file to parse. Use --lca_file=<file>\n";
}
if ( !$options{bam_file} && !$options{good_list} && !$options{bad_list} ) {
    die "ERROR: Must give reads to pull from the LCA file. Use either: --bam_file, --good_list, or --bad_list\n";
}

my %mapped_reads;

if ( $options{bam_file} ) {
    open( my $in, "-|", "samtools view -F0xC $options{bam_file}" ) or die "Unable to run samtools view on $options{bam_file}\n";
    print STDERR "Making a list of the M_M reads from the bam file . . .\n";
    while (<$in>) {
        my @fields = split( /\t/, $_ );
        $fields[0] =~ s/_(1|2)$//;
        $mapped_reads{ $fields[0] }++;
    }
    print STDERR "Finished reading in M_M reads from the bam.\n";
}

if ( $options{good_list} || $options{bad_list} ) {
    my $in;
    if ( $options{good_list} ) {
        open( $in, "<", $options{good_list} ) or die "Unable to open input good_list: $options{good_list}\n";
    }
    if ( $options{bad_list} ) {
        open( $in, "<", $options{bad_list} ) or die "Unable to open input bad_list: $options{bad_list}\n";
    }
    while (<$in>) {
        chomp;
        $_ =~ s/(_[12]{1})$//;
        $_ =~ s/(\s+)$//g;
        $mapped_reads{$_}++;
    }
    close $in;
    print STDERR "Finished reading in list of reads to parse.\n";
}

open( LCA, "<", $options{lca_file} ) or die "ERROR: Couldn't open the lca file: $options{lca_file}\n";
print STDERR "Starting to read through the LCA file, parsing for the desired reads . . .\n";
while (<LCA>) {
    chomp;
    my @line = split( /\t/, $_ );
    my $read = $line[0];
    if ( $read =~ /(_[12]{1})$/ ) { $read =~ s/(_[12]{1})$//; }
    if ( $options{bam_file} || $options{good_list} ) {
        if ( $mapped_reads{$read} ) {
            print STDOUT "$_\n";
        }
    }
    elsif ( $options{bad_list} ) {
        next if ( $mapped_reads{$read} );
        print STDOUT "$_\n";
    }
}
print STDERR "Finished parsing the LCA file.\n";

sub help {
    die "\nHELP: This script will parse an LCA file for the desired reads. Output is STDOUT.
   --lca_file=\tLCA file to pull the subset of lcas from.
   --bam_file=\tA bam file w/ \"good\" reads to pull from the lca file. Will only take M_M reads from the .bam
   --good_list=\tA list of reads that you want from the lca file.
   --bad_list=\tA list of reads to remove from the lca file\n\n";
}
