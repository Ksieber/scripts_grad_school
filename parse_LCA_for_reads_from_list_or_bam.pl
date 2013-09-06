#!/usr/local/bin/perl/ -w
use strict;

use Getopt::Long qw (:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(\%options,
               'lca_file=s',
               'bam_file=s',
               'good_list=s',
               'bad_list=s',
               'help'
      );

if($options{help}){die "HELP: This script will parse an LCA file for the desired reads. Output is STDOUT.
   --lca_file=\tLCA file to pull the subset of lcas from.
   --bam_file=\tA bam file w/ \"good\" reads to pull from the lca file. Will only take M_M reads from the .bam
   --good_list=\tA list of reads that you want from the lca file.
   --bad_list=\tA list of reads to remove from the lca file\n";
}

if(!$options{lca_file}){
   die "ERROR: You must give an LCA file to parse. Use --lca_file=<file>\n";
}
if(!$options{bam_file} && !$options{good_list} && !$options{bad_list}){
   die "ERROR: Must give reads to pull from the LCA file. Use either: --bam_file, --good_list, or --bad_list\n";
}

my %mapped_reads;   ## This is a terrible name for the hash. It is used to keep track of the reads to parse on. Ie good or bad reads depending on what is given.


if($options{bam_file}){
   open(my $in, "-|", "samtools view -F0xC $options{bam_file}") or die "Unable to run samtools view on $options{bam_file}\n";
   print STDERR "Making a list of the M_M reads from the bam file . . .\n";
   while(<$in>){
      my @fields = split(/\t/,$_);
      $fields[0]=~s/_(1|2)$//;   
      $mapped_reads{$fields[0]}++;
      #print "$fields[0]\t$mapped_reads{$fields[0]}\n";  ## Working
   }  
   print STDERR "Finished reading in M_M reads from the bam.\n";
}



if($options{good_list} || $options{bad_list}){
   my $in;
   if($options{good_list}){
      open($in, "<", $options{good_list}) or die "Unable to open input good_list: $options{good_list}\n";
   }
   if($options{bad_list}){
      open($in, "<", $options{bad_list}) or die "Unable to open input bad_list: $options{bad_list}\n";
   }
   while(<$in>){
      chomp;
      $_=~s/_(1|2)$//;  
      $mapped_reads{$_}++;
   }
   print STDERR "Finished reading in list of reads to parse.\n";
}



open(LCA, "<", $options{lca_file}) or die "ERROR: Couldn't open the lca file: $options{lca_file}\n";
print STDERR "Starting to read through the LCA file, parsing for the desired reads . . .\n";
while(<LCA>){
   chomp;
   my ($read,$lca)=split(/\t/,$_);
   $read=~s/_(1|2)$//;
#print STDERR "READ:$read\tLCA:$lca\n";  ## Working
   if($options{bam_file} || $options{good_list}){
      if($mapped_reads{$read}){
          print STDOUT "$_\n";
      } 
   }
  if($options{bad_list}){
    next if($mapped_reads{$read});
    print STDOUT "$_\n"; 
  }
}
print STDERR "Finished parsing the LCA file.\n";
   
