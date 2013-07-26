#!/usr/local/bin/perl -w
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
      'input=s',
      'dfam_acc=s',
      'output_dir=s',
      'help|h',
      );

if($options{help}){
   die "Help: This script will take the *.dfam_overlap output and convert the dfam #'s into the TE names.
      --input=\tThe name of the dfam_overlap output file.
      --dfam_acc=\tThe name of the file with the Dfam #'s and TE names.
      --output_dir=\tDirectory for output.
      --help\n";
}

if(!$options{input} || !$options{dfam_acc}){
   die "ERROR: Must give an input file and dfam_acc file. Use --input & --dfam_acc\n";
}

my %acc;

open(ACC, "<", "$options{dfam_acc}") or die "ERROR: Unable to open the dfam_acc file: $options{dfam_acc} because: $!\n";
while(<ACC>){
   chomp; 
   my @f=split;
   $acc{$f[0]}=$f[1];
}
close ACC or die "ERROR: Unable to close the dfam_acc file: $options{dfam_acc} because: $!\n";
my $output;
if($options{input}=~/(\w+?)\.dfam_overlap$/){
   $output=$1;
} else {
   $output="dfam_hits_with_names"
}

open(IN, "<", "$options{input}") or die "ERROR: Unable to open input: $options{input} because: $!\n";
open(OUT, ">", "$options{output_dir}/$output\.dfam_overlap_names") or die "ERROR: Unable to open output: $options{output_dir}/$output\.dfam_overlap_names\n";
while(<IN>){
   chomp; 
   my @q=split(/\t/);
   print OUT "$q[0]";
   shift @q;
   if($q[0] =~ /none/){
      print OUT "\tnone\n";
      next;
   }
   foreach my $hits (@q){
      if(!$acc{$hits}){
         print OUT "\t$hits"
      } else {
         print OUT "\t$acc{$hits}";
      }  
   }
  print OUT "\n";
}
close IN or die "ERROR: Unable to close input: $options{input} because: $!\n";
close OUT or die "ERROR: Unable to close output: $options{output_dir}/$output because: $!\n";








