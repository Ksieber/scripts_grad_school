#!/usr/local/bin/perl -w

##NOT WORKING YET
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options;
my $results = GetOptions (\%options,
		'input=s',
		'help|h',
		);

if ($options{help}){die "\nHELP: This script will\n";}
if (!$options{input} && !$ARGV[0]){die "Must give an input file with --input or by command line ARGV[0].  Try again.\n";}


my $input = $options{input} ? $options{input} : $ARGV[0];
my @rows;
my @transposed;

open(IN,"<","$input") or die "Can not open the input: $input because: $!\n";
while(<IN>){
   chomp; 
   my @row=split(/\s+/,$_);
   push(@rows,[@row]);
}

for my $row (@rows) {
   for my $column (0 .. $#{$row}) {
       push(@{$transposed[$column]}, $row->[$column]);
   }
}
  
for my $new_row (@transposed) {
   for my $new_col (@{$new_row}) {
       print $new_col, "\t";
   }
   print "\n";
}

