#!/usr/local/bin/perl5 -w
use strict;
use File::Basename;
if(@ARGV!=1){die "Must pass ARGV an input fasta file to substitute the id's for.\n";}
my ($file,$path,$suf) = fileparse($ARGV[0],".fa");
open(OUT,">","$path/$file\_sub-id.fa");
open(ID,">","id_list.txt");
open(IN,"<","$ARGV[0]");
my $last_read_seen=0;
my %foo;
my $a;

while(<IN>){
   chomp; 
   if($last_read_seen==1){
      $last_read_seen=0;
      next;
   }
   if($_!~/^\>/){
      print OUT "$_\n";
      next;
   }
   $_=~/^>(.*)/;
#print STDERR "$1\n";
   my $id=$1;
#print STDERR "$a\n";
   $id=~s/ /_/g;
   $id=~s/complete_genome//g;
   $id=~s/\,//g;
   $id=~s/\.//g;
   $id=~s/_$//;
   if($foo{$id}){
#print STDERR "THIS READ IS A DUP!: $id\n";
      $last_read_seen=1;
      next;
   }
   $a++;
   $foo{$id}="XYZ$a";
#print STDERR "$foo{$id}\t$id\n";
   print OUT "\>$foo{$id}\n";
   print ID "$foo{$id}\t$id\n";
}
close IN;
close ID;
close OUT;
