#!/usr/local/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
use warnings;
use strict;
use LGTSeek;
use File::Basename;
use run_cmd;
use POSIX;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
      'input=s',
      'threads=s',
      'output_prefix=s',
      'output_dir=s',
      'Qsub=s',
      'sort_mem=s',
      'help',
);

if($options{help}){
   die "Help: This script will take a bam to sort and index it.
      --input=           	Position sorted bam to index. May take one input from --input or ARGV
      --output_prefix=  	\$Prefix.bai
      --output_dir=     	/Directory/to/put/output
      --Qsub=			<0|1> [0] 1= Qsub to grid.
      --help\n";
}

if(!$options{input} && !$ARGV[0]){ die "Must pass an input bam to sort, use --input=<BAM> or pass \$ARGV[0]\n"; }
my $input = $options{input} ? $options{input} : $ARGV[0];
chomp($input);
my $lgtseek = LGTSeek->new2(\%options);
my ($bam,$path,$suf)=fileparse($input,'.bam');  ## @{$lgtseek->{bam_suffix_list}}
$options{output_prefix} = $options{output_prefix} ? $options{output_prefix} : $bam;
$options{output_dir} = $options{output_dir} ? $options{output_dir} : $path;
my $out = "$options{output_dir}\/$options{output_prefix}";
my $Qsub = defined $options{Qsub} ? $options{Qsub} : 0;

my $cmd= "samtools index $input $out\.bai";

if($Qsub==1){
      Qsub({
            cmd => $cmd,
            sub_name => "indexbam"
            })
} else {
      print STDERR "+++ Indexing bam +++\n";
      run_cmd($cmd);
      print STDERR "+++ Completed indexing: $input\n"; 
}


