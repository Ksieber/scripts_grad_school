#!/usr/local/bin/perl
use strict;
use warnings;
use setup_input;
use setup_output;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
our $results = GetOptions (\%options,
		'input=s',
		'input_list=s',
		'output_dir=s',
      'subdirs=s',
		'password=s',
      'test_print',
		'help',
);

if($options{help}){die
	"Help: This script will decrypt SRA data (.ncbi_enc)
	--input=
	--input_list=
	--output_dir=		[cwd]
	--password=
	--help\n";
}

if(!$options{input} && !$options{input_list}){die
	"Error: Must give an input. use --input= or --input_list=.\n";
}
if(!$options{password}){ die
   "Error: Must give the password to decrypt the files with --password=\n";
}

my $input=setup_input();           ## Returns array(ref) of input(s)
my $out_dir=setup_output($input);   ## Returns hash(ref) of output dirs Key = input, value = output_dir
foreach my $file (@$input){
   if($options{test_print}){
      print STDERR "/local/projects-t3/HLGT/TCGA/decrypt/decrypt.bin $file -password $options{password} -out-dir $out_dir->{$file}->{dir}\n";
   } else { 
      `/local/projects-t3/HLGT/TCGA/decrypt/decrypt.bin $file -password $options{password} -out-dir $out_dir->{$file}-{dir}`;
   }
}

