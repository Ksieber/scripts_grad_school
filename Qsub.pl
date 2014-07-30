#!/usr/bin/perl
use warnings;
no warnings 'uninitialized';
use strict;
use run_cmd;
use lib ('/home/ksieber/scripts/','/local/projects-t3/HLGT/scripts/lgtseek/') or die "Error: Unrecognized command line option. Please try again.\n";
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
	'cmd=s', 	
	'log=s',
	'name=s',
	'threads=s',
	'mem=s',
	'project=s',
	'cwd=s',
	'wd=s',
	'help'
	);
if($options{help}){die 
	"This script will Qsub a command. The following options can be passed.
	--cmd 		= command/shell script to be qsub'ed [STDIN].
	--log 		= file to log in.
	--name		= Job name [STDIN].
	--threads	= # cpu threads to use [1].
	--mem		= min amount of RAM to use.
	--project	= grid project to use [jdhotopp-lab].
	--cwd		= <0|1> [0] 1= Use current working directory. 
	--wd 		= Directory for grid to work from [~/].\n";
}

while(<>){
	chomp;
	my $cmd = defined $options{cmd} ? "$options{cmd}" : "$_";
	my $report = Qsub2({ 
		cmd => $cmd, 
		log => $options{log}, 
		name => $options{name}, 
		threads => $options{threads}, 
		mem => $options{mem},
		project => $options{project},
		cwd => $options{cwd},
		wd => $options{wd}
		});
}

