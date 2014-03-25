#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
=head1 NAME

blast2lca.pl

=head1 SYNOPSIS

Calculate the LCA from hits in a blast-m8 report.

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
		'fasta|f=s',
		'bam|b=s',
		'db=s',
		'formatdb=i',
		'Qsub|Q=i',
		'sub_name=s',
		'print_hostname|ph=i',
		'output_dir|o=s',
		'output_prefix=s',
		'subdirs=i',
		'conf_file=s',
		'verbose=i',
		'help|h',
		) or die "Error: Unrecognized command line option. Please try again.\n";

use lib ("/local/projects-t3/HLGT/scripts/lgtseek/lib/","/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/");
use print_call;
print_hostname(\%options);						## This is useful for trouble shooting grid nodes that might be missing modules for LGTSeek etc. 

use LGTSeek;
use run_cmd;
use setup_input;
use Carp;
$Carp::MaxArgLen = 0;
use File::Basename;

if($options{help}){&help;}		## @ end of script
if(!$options{fasta} && !$options{bam}){confess "Error: Please give an input with --fasta=<FASTA> or --bam=<BAM>.\n";}

## Print the script call
print_call(\%options);
## Initialize LGTSeek.pm
my $lgtseek = LGTSeek->new2(\%options);

## Setup output directory
my $input = $options{fasta} ? $options{fasta} : $options{bam};
my ($name,$path,$suf) = fileparse($input,$lgtseek->{suffix_regex});
chomp $name;
if(!$lgtseek->{output_dir}){$lgtseek->{output_dir}="$path/lgtseq/";}
if($lgtseek->{subdirs}){
	$lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}"); 
	$lgtseek->{output_dir} = "$options{output_dir}/"."$name/"; 
	$lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}");
}
if($options{formatdb}){
	run_cmd("formatdb -p F -i $options{db}");
}
# Run Blast
my $blast;
if($options{fasta}){
	$blast = $lgtseek->bestBlast2({
		fasta => $input,
		db => $options{db},
		output_dir => $lgtseek->{output_dir},
		});
} else {
	$blast = $lgtseek->bestBlast2({
		bam => $input,
		db => $options{db},
		output_dir => $lgtseek->{output_dir},
		});
}

print_complete(\%options);

sub help {
	die "Help: This script will take a blast-m8 report and calculate the LCA for each unique ID.
	----------------------------------------------------------------------------------------------------
		--fasta|f=			<FASTA> File to blastn against the db. 
		--bam|b=			 <BAM>  File to blastn against the db. 
	----------------------------------------------------------------------------------------------------
		--db|d=				<FastaDB> Fasta reference file to blastn against. 
		--formatdb=			<0|1> [0] 1= Build the blast database to blast again. 
	----------------------------------------------------------------------------------------------------
		--output_dir|o=			Directory for all output. Will be created if it doesn't exist.
		  --output_prefix=		Name of prefix for the output. 
		  --subdirs=			<0|1> [0] 1 = Create a new subdirectory for each input in the output directory.
	----------------------------------------------------------------------------------------------------
		--Qsub|Q=			<0|1> [0] 1 = Submit the job to the grid.
		  --sub_name=  			Name of the job to submit.
		  --project=			[jdhotopp-lab] Grid project to use.
	----------------------------------------------------------------------------------------------------	
		--conf_file=			[/home/USER/.lgtseek.conf]
	----------------------------------------------------------------------------------------------------	
		--verbose|V			<0|1> [1] 1= Verbose reporting of progress. 0 =Turns off reports.
		--print_hostname|ph= 		<0|1> [1] Print hostname. Defaults to 1 if verbose=1. 
		--help
	----------------------------------------------------------------------------------------------------\n";
}
















