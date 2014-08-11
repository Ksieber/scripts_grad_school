#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
=head1 NAME

blast2lca.pl

=head1 SYNOPSIS

Append Taxonomy to blast-m8 report.

=head1 DESCRIPTION


=head1 AUTHORS - Karsten Sieber

e-mail: Karsten.sieber@gmail.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with "_"

=cut
use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
		'input|i=s', # Comma separated list of files
		'input_list|I=s',
		'Qsub|Q=i',
		'sub_name=s',
		'project=s',
		'print_hostname|ph=i',
		'output_dir|o=s',
		'output_prefix=s',
		'subdirs=i',
		'taxon_host=s',
		'taxon_dir=s',
		'taxon_idx_dir=s',
		'conf_file=s',
		'verbose=i',
		'help|h',
		) or die "Error: Unrecognized command line option. Please try again.\n";

use run_cmd;
if($options{Qsub}){ &Qsub_script(\%options) };
use print_call;
print_hostname(\%options);									## This is useful for trouble shooting grid nodes that might be missing modules for LGTSeek etc. 

use lib ('/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/','/local/projects-t3/HLGT/scripts/lgtseek/lib/');
use LGTSeek;
use setup_input;
use Carp;
$Carp::MaxArgLen = 0;
use File::Basename;

if($options{help}){&help;}		## @ end of script
if(!$options{input} && !$options{input_list}){confess "Error: Please give an input.bam with --input=<FILE> or --input_list=<LIST>. Try again or use --help_full.\n";}

## Print the script call
print_call(\%options);
## Initialize LGTSeek.pm
my $lgtseek = LGTSeek->new2(\%options);
## Setup array ref of inputs
my $inputs = setup_input(\%options);
# Setup a few default values
my $evalue_cutoff = $options{evalue_cutoff} ? $options{evalue_cutoff} : "1";
my $bho = $options{best_hits_only} ? $options{best_hits_only} : "0";
my $last_id;
# my $combine_PE_lca= $options{combine_PE_lca} ? $options{combine_PE_lca} : "0";	## NOT implemented yet

# Intialize gi2taxon db
my $gi2tax = $lgtseek->getGiTaxon();

#
foreach my $input (@$inputs){
	my ($name,$path,$suf) = fileparse($input,$lgtseek->{suffix_regex});
	chomp $name;
	open(IN, "<$input") or confess "Error: Unable to open input: $input .\n";

	## Setup output directory

	my $out_dir = $options{output_dir} ? $options{output_dir} : "$path\/blast2tax/";
	if($lgtseek->{subdirs}){
		$lgtseek->_run_cmd("mkdir -p $out_dir"); 
		$out_dir = "$out_dir/"."$name/"; 
		$lgtseek->_run_cmd("mkdir -p $out_dir");
	}
	my $out_pref = $options{output_prefix} ? $options{output_prefix} : $name;
	my $output;
	if($options{output_dir} || $options{output_prefix}){
		$lgtseek->_run_cmd("mkdir -p $out_dir"); 
		open($output, ">", "$out_dir\/$out_pref\_with-taxonomy.txt") or confess "Error: Can't open output: $out_dir\/$out_pref\_with-taxonomy.txt because: $!\n";
	} else { 
		$output = *STDOUT;
	}

	if($lgtseek->{subdirs}){
		$lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}"); 
		$lgtseek->{output_dir} = "$options{output_dir}/"."$name/"; 
		$lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}");
	}

	while(<IN>){
		chomp; 
		my @fields = split;
		my $taxon = $gi2tax->getTaxon($fields[1]);
		print $output "$_ $taxon->{lineage}\n";	
	}
	close $output;
}

print_complete(\%options);

sub help {
	die "Help: This script will take a blast-m8 report and calculate the LCA for each unique ID.
	----------------------------------------------------------------------------------------------------
		--input|i=			<Input> Accepts blast -m8 reports.
		--input_list=			List of blast -m8 reports, 1 file per line.
	----------------------------------------------------------------------------------------------------
		*** Output will go to STDOUT unless given an --output_dir or --output_prefix ***
		--output_dir|o=			Directory for all output. Will be created if it doesn't exist.
		  --output_prefix=		Name of prefix for the output. 
		  --subdirs=			<0|1> [0] 1 = Create a new subdirectory for each input in the output directory.
	----------------------------------------------------------------------------------------------------
		--Qsub|Q=			<0|1> [0] 1 = Submit the job to the grid.
		  --sub_name=  			[STDIN] Name of the job to submit.
		  --project=			[jdhotopp-lab] Grid project to use. 
	----------------------------------------------------------------------------------------------------	
		--taxon_host=			[mongotest1-lx.igs.umaryland.edu:10001]
		--taxon_dir=			[/local/db/repository/ncbi/blast/20120414_001321/taxonomy/]
		--taxon_idx_dir=		[/local/projects-t3/HLGT/idx_dir/20120414]
	----------------------------------------------------------------------------------------------------	
		--print_hostname|ph= 		<0|1> [0] Print hostname. Defaults to 1 if verbose=1. 
		--conf_file=			[/home/USER/.lgtseek.conf]
		--help
	----------------------------------------------------------------------------------------------------\n";
}
















