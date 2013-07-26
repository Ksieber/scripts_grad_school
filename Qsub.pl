#!/usr/local/bin/perl -w
use strict;

while(<>){
   chomp; 
   my $report=`qsub -V -P jdhotopp-lab $_`;
   chomp($report);
   print STDERR "$report\n";
}
