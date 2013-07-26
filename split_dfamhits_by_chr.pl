#!/usr/local/bin/perl -w
#use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
                          'input=s',
                          'help',
                         );

if($options{help}){die "HELP: This script will split up the Dfam_hits file into each chromosome. Use --input.\n";}

my %hash;
my $reference;

open(IN, "<", $options{input}) or die "Unable to open input file: $options{input} because :$!\n";
while(<IN>){
    chomp;
    my @f=split(/\t/,$_);
    my @hit = split(/\./,$f[15]);
    if($hit[0] =~ /^GL(\d+)/){
        $reference = "chrUn_gl$1";
    } elsif ($hit[0] =~ /^(\d+|X|Y)/){
        $reference = "chr$1";
    } 
#   print STDERR "$reference\n";
    if($hash{$reference}){x
        print $reference "$_\n"; 
    }
    if(!$hash{$reference}){
        open($reference, ">", "$reference\.Dfam_hits") or die "Unable to open output: $reference\.Dfam_hits because: $!\n";
        print $reference "$_\n";
        $hash{$reference}++;    
    }
}
