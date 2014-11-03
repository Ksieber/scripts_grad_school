#!/usr/bin/perl -w
use strict;
use run_cmd;
use lib ('/home/ksieber/scripts/','/local/projects-t3/HLGT/scripts/lgtseek/');


if($ARGV[0]){die "Error: Pipe this script a qsub command or shell script through STDIN\n";}

while(<>){
	chomp;
   	my $report=Qsub("$_");
}
