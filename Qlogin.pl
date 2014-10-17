#!/usr/bin/perl
use warnings;
no warnings 'uninitialized';
use strict;
use run_cmd;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'threads=s', 'sub_mem=s', 'project=s', 'help' ) or die "Error: Unrecognized command line option. Please try again.\n";

my $project = defined $options{project} ? " -P $options{project}"         : " -P jdhotopp-lab";
my $threads = defined $options{threads} ? " -pe thread $options{threads}" : undef;
my $sub_mem = defined $options{sub_mem} ? " -l mf=$options{sub_mem}"      : " -l mf=8G";
my $cmd     = "qlogin$project$threads$sub_mem -q interactive.q";

system($cmd);
