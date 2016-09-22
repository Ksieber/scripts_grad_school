#!/usr/bin/perl
use warnings;
use strict;
use run_cmd;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'mem|m=s', 'Picard_jar=s', 'argv_string=s', 'help|?', 'Qsub|q=i', 'sub_mem=s', );

if ( $options{help} ) {
    die "This script will launch Picard.jar
    Arguements for Picard are passed through ARGV or --argv_string=\"string\". ex: picard.pl SamToFastq I=input F=out.fq
    To get help information for each component of Picard use plain \"help\" after the component name. ex: picard.pl SamToFastq help
	--mem= 			[2g] Amount of mem for java to use.
	--Qsub=			<0|1> [0] 1= Qsub.
	  --sub_mem=	[3G]
	  --sub_mail=	[0]
	  --wd=			
	--Picard_jar=	[/home/ksieber/lib/picard/dist/picard.jar]
	--java=			[java]\n";
}

$options{argv_string} = defined $options{argv_string} ? $options{argv_string} : join( " ", @ARGV );

my $mem         = defined $options{mem}         ? $options{mem}         : "2g";
my $Picard_jar  = defined $options{Picard_jar}  ? $options{Picard_jar}  : "/home/ksieber/lib/picard/dist/picard.jar";
my $java        = defined $options{java}        ? $options{java}        : "java";
my $argv_string = defined $options{argv_string} ? $options{argv_string} : join( " ", @ARGV );
my $sub_mem     = defined $options{sub_mem}     ? $options{sub_mem}     : "3G";

my $cmd = "$java -Xmx$mem -XX:-UseGCOverheadLimit -XX:-UseParallelGC -jar $Picard_jar $argv_string";

if ( $options{Qsub} ) { Qsub({ cmd => $cmd, sub_mem => $sub_mem }); }
else                  { `$cmd`; }

