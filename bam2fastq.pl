#!/usr/local/bin/perl 
use warnings;
use strict;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
      'input=s',
      'fasta=s',
      'no_nums=s',
      'interleaved=s',
      'help|h',
      );
if($options{help}){die
   "This script will take a bam and output the fastq. 
    --input=		input bam or use ARGV[0].
    --fasta= 		<0|1> [0] 0=Fastq. 1=fasta
    --interleaved= 	<0|1> [0] 0=Seperate files for output(_1/_2). 1=Interleaved pairs.
    --no_nums=  	<0|1> [0] 0=Add _1 and _2 to the fastq. 1=Don't add numbers.\n";
}
## Open Input
if(!$ARGV[0] && !$options{input}){die
   "Must give an input bam with ARGV[0] or --bam=<input>.\n";
}
my $in = $options{input} ? $options{input} : $ARGV[0];
open(IN,"-|","samtools view $in") || die "Unable to open input $in.\n";
my ($fn,$path,$suf)=fileparse($in,'.bam');

## Set Defaults
if(!$options{no_nums}){$options{no_nums}=0;}
if(!$options{fasta}){$options{fasta}=0;}
if(!$options{interleaved}){$options{interleaved}=0;}
if($options{no_nums}==1){$options{interleaved}=0;} 						## If _1/_2 are NOT added, MUST prints reads to different files

## Setup output
my $out1;
my $out2;
if ($options{interleaved}==0) {
	if ($options{fasta}==0) {
		open($out1,">","$path/$fn\_1.fq") or die "Can't open output file: $path/$fn\_1.fq because: $!\n";
		open($out2,">","$path/$fn\_2.fq") or die "Can't open output file: $path/$fn\_2.fq because: $!\n";
	} elsif ($options{fasta}==1) {
		open($out1,">","$path/$fn\_1.fa") or die "Can't open output file: $path/$fn\_1.fa because: $!\n";
		open($out2,">","$path/$fn\_2.fa") or die "Can't open output file: $path/$fn\_2.fa because: $!\n";
	}
} elsif ($options{interleaved}==1) {
	if ($options{fasta}==0) {
		open($out1,">","$path/$fn\.fq") or die "Can't open output file: $path/$fn\.fq because: $!\n";
		$out2=$out1;
	} elsif ($options{fasta}==1) {
		open($out1,">","$path/$fn\.fa") or die "Can't open output file: $path/$fn\.fa because: $!\n";
		$out2=$out1;
	}
}

## Process .bam conversion to .fastq
my %foo;
while(<IN>){
	chomp;
	next if $_=~/^@/;
	my @f=split;
   	$foo{$f[0]}++;
   	if($foo{$f[0]} == 2){
		if($options{no_nums}==0){$f[0] =~ s/$f[0]/$f[0]\_2/;}
        _print($f[0],$f[9],$f[10],$out2);
    	next;
    } 
    if($foo{$f[0]} == 1){
    	if($options{no_nums}==0){$f[0] =~ s/$f[0]/$f[0]\_1/;}
      	_print($f[0],$f[9],$f[10],$out1);
   	}	
}

sub _print {
	my ($read,$seq,$qual,$fh) = @_;
	if($options{fasta}==0){print $fh "\@$read\n$seq\n\+\n$qual\n";}
	if($options{fasta}==1){print $fh "\>$read\n$seq\n";}
}

  
  
  
  
  
  
