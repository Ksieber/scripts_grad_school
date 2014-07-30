#!/usr/local/bin/perl -w
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options;
my $results = GetOptions (\%options,
			  'ref_list=s',
			  'output_dir=s',
			  'refloop_list',
			  'help|h',
			 );
if ($options{help}){die "\nThis script will take a list of references and bwa index them.  
\t ref_list=list of reference files to index
\t output_dir= directory where to save the indexing\n\n";}

my @index_list;

open (REFLIST, "<", $options{ref_list}) || die "ERROR: Couldn't open $options{ref_list} becaues $!\n";
while (<REFLIST>){
  push(@index_list, $_);
}

if ($options{refloop_list}){
  open (REFLOOP, ">", "/local/projects/HLGT/ksieber_dir/BWA_indexed.list") || die "ERROR: Couldn't open output list because $!\n";
}

foreach my $file (@index_list){
  $file =~ m/(\w{2}_\d+)\.?\d?/;
  print STDERR "\nProcesssing $file";
  `bwa index -p $options{output_dir}$1 $file`;
  if ($options{refloop_list}){
    print REFLOOP "$1\n";
  }  
}



__END__
