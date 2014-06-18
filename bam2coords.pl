#!/usr/bin/perl

use strict;
use lib '/local/projects-t3/HLGT/scripts/lgtseek-master/lib';
use LGTSeek;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
                          'input_bam=s', # Comma separated list of files
                          'good_list=s',                        
                          'output_dir=s',
                          'output_prefix=s',
             			  'help|h',
);

if($options{help}){die "HELP: This script will take a bam and report the chromosome, cordinates and orientation of each mapped read. 
	--input_bam=s		<Input bam>. Required.
	--good_list=s		A list of reads to pull cords for. Reads list into memory.
	--output_dir=s		Directory for output. [input dir]
	--output_prefix=s	Prefix for the output file. [input bam]
	--bin_dir=s		Directory with LGTSeek.pm [/local/projects-t3/HLGT/scripts/lgtseek-master/lib/]
	\n";
}

if(!$options{input_bam}){die "Error: Please give an input.bam with --input_bam=<FILE>. Try again or use --help.\n";}

# Take care of the inputs
my $bin_dir = $options{bin_dir} ? $options{bin_dir} : '/local/projects-t3/HLGT/scripts/lgtseek-master/lib/';

my ($fn,$dir,$suf) = fileparse($options{input_bam},".bam");
#$dir =~s/\/$//;
my $out_prefix = $options{output_prefix} ? $options{output_prefix} : $fn;
my $out_dir = $options{output_dir} ? $options{output_dir} : $dir;
my $out = "$out_dir\/$out_prefix\.txt";

# Create an lgtseek object
my $lgtseek = LGTSeek->new({
    bin_dir => $bin_dir,
    output_dir => $options{output_dir},
    paired_end => 1,
});

my %good_list;

if($options{good_list}){
	open(LIST,"<","$options{good_list}") or die "Can't open the good_list: $options{good_list} because: $!\n";
	while(<LIST>){
		chomp;
		$good_list{$_}++;
	}
	close LIST;
}

open(OUT,">","$out") or die "Can't open output: $out because: $!\n";
print OUT "Reference_name\tPosition\tbp\tStrand\n";
open(IN,"samtools view -h $options{input_bam} |") or die "Can't open input: $options{input_bam} because: $!\n";
while(<IN>){
	chomp;
	my ($qname,$flag,$rname,$pos,$mapq,$cigar,$mrnm,$mpos,$tlen,$seq,$qual,$opt)=split;
	if($good_list{$qname}){
		my $strand;
		my $parseFlag = $lgtseek->_parseFlag($flag);
		next if($parseFlag->{'qunmapped'}==1);
		my $length=length($seq);
		if($parseFlag->{'qrev'}==1){
			$strand = '-';
		} else {
			$strand = '+';
		}
		print "$rname\t$pos\t$length\t$strand\n";
	}
}
close OUT;
close IN;
	


	
