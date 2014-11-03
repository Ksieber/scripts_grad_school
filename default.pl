#!/usr/local/bin/perl -w
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
		'input=s',
		'input_list=s',
		'output_dir=s',
		'password=s',
		'help',
      );
if($options{help}){die
	"Help: This script will decrypt TCGA data. 
	--input=
	--input_list=
	--output_dir=		[cwd]
	--password=
	--help\n";
}

if(!$options{input} && !$options{input_list}){die
	"Error: Must give an input. use --input= or --input_list=.\n";
}

