#!/usr/local/bin/perl 
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
      'bam=s',
      'fasta',
      'add_num',
      'help|h',
      );
if($options{help}){die
   "This script will take a bam and output the fastq to STDOUT. samtoosl bam2fq is faster?
    Use either ARGV[0]
    --bam= to give input bam.
    --add_num to add _1 and _2 to the fastq.
    --fasta for fasta output instead of fastq.\n";
}
if(!$ARGV[0] && !$options{bam}){die
   "Must give an input bam with ARGV[0] or --bam=<input>.\n";
}
my $in = $options{bam} ? $options{bam} : $ARGV[0];
open(IN,"-|","samtools view $in") || die "Unable to open input $in.\n";
my %foo;
while(<IN>){
   chomp;
   next if $_=~/^@/;
   my @f=split;
   if($options{add_num}){
      $foo{$f[0]}++;
      if($foo{$f[0]} == 2){
         $f[0] =~ s/$f[0]/$f[0]\_2/;
      } 
      if($foo{$f[0]} == 1){
         $f[0] =~ s/$f[0]/$f[0]\_1/;
      }
   }
  if($options{fasta}){
    print STDOUT "\>$f[0]\n$f[9]\n";
    next;
  } 
  print STDOUT "\@$f[0]\n$f[9]\n\+\n$f[10]\n";
}
