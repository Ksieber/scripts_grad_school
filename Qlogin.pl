#!/usr/bin/perl
use warnings;
no warnings 'uninitialized';
use strict;
use run_cmd;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'threads|t=s', 'sub_mem|m=s', 'project|p=s', 'help|?' ) or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) {
    die "This script will Qlogin into the SGE grid. 
	--threads|t=		< # > Number of threads to use. [1]
	--sub_mem|m=		< #G > Amount of needed memory. [8G]
	--project|p=		SGE project group. 		[jdhotopp-lab]
	--help|?		This wonderful helpful information.\n";
}

my $project = defined $options{project} ? " -P $options{project}"         : " -P jdhotopp-lab";
my $threads = defined $options{threads} ? " -pe thread $options{threads}" : undef;
my $sub_mem = defined $options{sub_mem} ? " -l mf=$options{sub_mem}"      : " -l mf=8G";
my $cmd     = "qlogin$project$threads$sub_mem -q interactive.q";

system($cmd);
