#!/usr/local/bin/perl -w
use strict;
use lib '/home/ksieber/scripts/';
use run_cmd;
use setup_input;
use setup_output;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
our $results = GetOptions (\%options,
   	      'input=s',
   	      'input_list',
	      'output_dir=s',
	      'cmd_log=s',
	      'Qsub=s',
	      'help',
);

if ($options{help}){die "\nHELP: This script will use picard tools to calculate the insert size for a sorted.bam.
\t--input=			REQUIRED. Must be a sorted.bam
\t--input_list=			List of sort.bams
\t--output_dir=       		Directory where to write the data if not the current working directory.
\t--cmd_log=			<0|1> [0] 1=print commands run to log file.
\t--Qsub=				<0|1> [0] 1=qsub each inputs Calculation\n";
}

if (!$options{input} && !$options{input_list}){die "ERROR: Must give a input.  Use: --input=some.sorted.bam or --input_list=<LIST of sorted bams>\n";}

## Setup input file(s)
my $input = setup_input();

## Setup outout directories
my $out_dir = setup_output($input);

## Setup Command logs
my $cmd_logs; 
if ($options{cmd_log}==1) {
	$cmd_logs=setup_logs($out_dir);
}


## Run Insert Calculation
foreach my $file (@$input){
	my $log;
	if($options{cmd_log}==1){
		$log=$cmd_logs->{$file};
	}
	my($fn,$path,$suf)=fileparse($file,('.srt.bam', '.bam', 'sort.bam'));
	my $out = "$out_dir->{$file}/$fn";
	my $cmd1="java -Xmx3g -jar /home/jdhotopp/bin/Picard/picard-tools-1.48/CollectInsertSizeMetrics.jar AS=true I=$file H=$out.std_insert.histogram O=$out.std_insert.metrics M=0 VALIDATION_STRINGENCY=SILENT DEVIATIONS=20";
	my $cmd2="java -Xmx3g -jar /home/jdhotopp/bin/Picard/picard-tools-1.48/CollectInsertSizeMetrics.jar AS=true I=$file H=$out.lrg_insert.histogram O=$out.lrg_insert.metrics M=0 VALIDATION_STRINGENCY=SILENT DEVIATIONS=1000000000000000000";
	if($options{Qsub}==1){
		Qsub($cmd1,$log);
		Qsub($cmd2,$log);
	} else {
		run_cmd($cmd1,$log);
		run_cmd($cmd2,$log);
	}
}

