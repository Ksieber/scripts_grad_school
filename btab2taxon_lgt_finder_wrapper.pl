#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use LGTBestBlast;

my %options;
my $results = GetOptions(\%options,
  'input=s',
  'output_dir=s',
  'output_prefix=s',
  'dbhost=s',
  'idx_dir=s',
  'taxon_dir=s',
  'min_length=s',
  'max_overlap=s',
  'help|h',
);

if($options{help}){die 
"Help: This script will take a blast m8 and add taxon info and search it for LGT with lgt_finder.
\t--input=\tFasta file m8 format. REQUIRED.
\t--output_dir=\tDirectory for output.
\t--output_prefix=\tPrefix for output.
\t--dbhost=
\t--idx_dir=
\t--taxon_dir=
\t--min_length=\tMin Length of blast hit for lgt_finder.
\t--max_overlap=\tMax overlap allowed in lgt_finder.\n";
}

my @suffix_list=('.txt');
my($file,$path,$suffx)=fileparse($options{input});
my $out_dir = $options{output_dir} ? $options{output_dir} : $path;
my $out_pref = $options{output_prefix} ? $options{output_prefix} : $file;
my $dbhost = $options{dbhost} ? $options{dbhost} : "mongotest1-lx.igs.umaryland.edu";
my $idx_dir = $options{idx_dir} ? $options{idx_dir} : "/local/projects-t3/HLGT/idx_dir/20120414/";
my $taxon_dir = $options{taxon_dir} ? $options{taxon_dir} : "/local/db/repository/ncbi/blast/20120414_001321/taxonomy";
my $min = $options{min_length} ? $options{min_length} : 1;
my $max = $options{max_overlap} ? $options{max_overlap} : 20;

`perl /home/ksieber/scripts/lgt_blastBesthit.pl --input=$options{input} --dbhost=$dbhost --taxon_dir=$taxon_dir --idx_dir=$idx_dir --overalloutput=$out_dir\/$out_pref\_besthits_taxon.txt`;
`perl /local/projects-t3/HLGT/scripts/lgt_finder_dr.pl --input=$out_dir\/$out_pref\_besthits_taxon.txt --output_prefix=$out_pref --min_length=$min --max_overlap=$max`;
