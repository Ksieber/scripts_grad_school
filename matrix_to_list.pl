#!/usr/local/bin/perl -w
use strict;
use Bio::Matrix::Generic;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

###This script will transform a matrix into a long list.

my %options;
my $results = GetOptions (\%options,
			  'input=s', #name of the file to transfrom
			  'output=s',
			  'help|h'
			 );
 
if ($options{help}){die "Help: This script will transform a file with a tab delimited matrix into a long list.
Use: --input= and --output=\n";}

##My Varialbes:
my $input=$options{input};
my $output=$options{output};
my $matrix = Bio::Matrix::Generix->new();
my $index;

open (INPUT, "<", $input) || die "ERROR: Couldn't open $input because $!\n";
open (OUTPUT, ">", $output) || die "ERROR: Couldn't open $output because $!\n";



while (<INPUT>){
  chomp;
  ##@line=split('\t',$_);
  $matrix->add_row($index, $_); 
  my @row = $matrix->get_row('A');
  print @row;
}


  
  
