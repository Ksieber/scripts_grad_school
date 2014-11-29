#!/usr/bin/perl
use warnings;
no warnings 'uninitialized';
use strict;
use run_cmd;
use lib ( '/home/ksieber/scripts/', '/local/projects-t3/HLGT/scripts/lgtseek/' );
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'cmd|c=s', 'log=s', 'sub_name|n=s', 'threads|t=s', 'sub_mem|mem=s', 'project|p=s', 'cwd=s', 'wd=s', 'sub_mail|ml=s', 'help|?' )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) {
    die "This script will Qsub a command. The following options can be passed.
    --cmd=          command/shell script to be qsub'ed [STDIN].
    --threads=      # cpu threads to use [1].
    --sub_name=         Job name [STDIN].
    --sub_mem=          min amount of RAM to use.
    --sub_mail=         1= email {\$user}\@som.umaryland.edu or the --sub_mail\@foo.bar     
    --log=          file to log in.
    --project=      grid project to use [jdhotopp-lab].
    --cwd=          <0|1> [0] 1= Use current working directory. 
    --wd=           Directory for grid to work from [~/].\n";
}

my $stdin;
if ( !$options{cmd} ) { chomp( $stdin = <> ); }
my $cmd = defined $options{cmd} ? "$options{cmd}" : "$stdin";
if ( $cmd !~ /\w+/ ) { die "Error: Unable to determine the cmd to qsub.\n"; }

my $report = Qsub(
    {   cmd      => $cmd,
        log      => $options{log},
        sub_name => $options{sub_name},
        sub_mem  => $options{sub_mem},
        sub_mail => $options{sub_mail},
        threads  => $options{threads},
        project  => $options{project},
        cwd      => $options{cwd},
        wd       => $options{wd}
    }
);

