#!/usr/local/bin/perl -w
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
      'blast=s',
      'output=s',
      'help',
      );

if($options{help}){die 
   "HELP:
      This script will take a blast-m8 file and parse the best hits (based on e-value). 
      --blast= Name of the input blast file. 
      --output= Name of the output file.\n";
}

if(!$options{blast} || !$options{output}){die
   "ERROR: Must give an input and output file. Please use --blast and --output\n";
}

open(IN, "<", $options{blast}) or die "Unable to open the blast file for parsing: $!\n";
open(OUT, ">", $options{output}) or die "Unable to open the output file: $!\n";

my %hash;

while(<IN>){
   chomp;
   my @f=split(/\t/,$_);
   my $read = $f[0];
   my $eval = $f[10];
#print STDERR "READ:$read\tEVAL$eval\n";
   if($hash{$read}){
#print STDERR "Hash exists ...\n";
#print STDERR "Checking if $eval is equal to $hash{$read}\n";
      if($eval eq $hash{$read}){
         #print STDERR "$eval equals $hash{$read}\n";
         print OUT "$_\n"; 
      }
   } else {
      if(!$hash{$read}){
         $hash{$read}=$eval;
         print OUT "$_\n";
      }
   }  
}

close IN;
close OUT;
