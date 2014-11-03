#!/usr/local/bin/perl5 -w
use strict;
use File::Basename;
if(@ARGV!=1){die "Must pass ARGV an input fasta file to substitute the id's for.\n";}
my ($file,$path,$suf) = fileparse($ARGV[0],".txt");
$path=~s/\/$//;
open(OUT,">","$path/$file\_sub-id.txt");
open(ID,"<","id_list.txt");
open(IN,"<","$ARGV[0]");
my %foo;

while(<ID>){
   chomp; 
   my @f=split;
   $foo{$f[0]}=$f[1];
}
close ID;

while(<IN>){
   $_=~s/(XYZ\d+)/$foo{$1}/g;
   print OUT;
} 
close IN;
close OUT;
