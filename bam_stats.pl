#!/usr/bin/perl

use strict;
use warnings;
use Time::SoFar;
use File::Basename;
use LGTSeek;
use setup_input;
use lib '/local/projects-t3/HLGT/scripts/lgtseek/lib/';
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
our $results = GetOptions (\%options,
	'input=s',
	'input_list=s',
	'output_dir=s',
	'subdirs=s',
	'append=s',
	'help'
);

if($options{help}){ die
	"This script will calculate the number of reads in a bam: Total, MM, MU, UU, SC
	--input=		<BAM>
	--input_list=	<list of BAMS>
	--output_dir=	</path/for/output/>
	--append=	</path/file/append.txt> Add output stats to this file.
	--help\n";
}

if(!$options{input} && !$options{input_list}){die "Must give an --input=<BAM> or --input_list=<list_of_bams>\n";}
my $append = $options{append} ? $options{append} : 0;

my $lgtseek = LGTSeek->new2();

my $input = setup_input();

foreach my $bam (@$input){
	if($lgtseek->empty_chk({input => $bam}) ==0){print STDERR "This is an empty bam! $bam\n"; next;}
	my $start = Time::SoFar::runinterval();
	open(my $in,"-|","samtools view $bam") or die "Can't open input: $bam\n";
	## 
	my @read_types = ('Total','MM','MU','UU','SC');
	my %counts;
	map { $counts{$_}=0; } @read_types;
	while(<$in>){
		chomp;
		my ($flag,$cigar) = (split /\t/, $_)[1,5];
		my $converted_flag=$lgtseek->_parseFlag($flag);
		if($cigar =~ /(\d+)M(\d+)S/ && $2 >= 24){$counts{SC}++;}
		if($cigar =~ /(\d+)S(\d+)M/ && $1 >= 24){$counts{SC}++;}
		if(!$converted_flag->{'qunmapped'} && !$converted_flag->{'munmapped'}){$counts{MM}++;}
		if(!$converted_flag->{'qunmapped'} && $converted_flag->{'munmapped'}){$counts{MU}++;}
		if($converted_flag->{'qunmapped'} && !$converted_flag->{'munmapped'}){$counts{MU}++;}
		if($converted_flag->{'qunmapped'} && $converted_flag->{'munmapped'}){$counts{UU}++;}
		$counts{Total}++;
	}
	## 
	my $OUT;
	my $print_header;
	if($options{output_dir}){
		my ($fn,$path,$suf)=fileparse($bam,".bam");
		my $out = $options{output_dir} ? "$options{output_dir}/$fn.stats" : "$path$fn.stats";
		open($OUT,"> $out") or die "Error: Can not open: $out\n";
		$print_header=1;
	} elsif ($options{append}){
		if(-e $options{append}){$print_header=0;} 
		else {$print_header=1;}
		open($OUT,">> $options{append}") or die "Error: Can not open: $append\n";
	} else {
		$OUT = *STDOUT;
		$print_header=1;
	}
	##
	if($print_header==1){print $OUT "BAM\t"; print $OUT join("\t",@read_types)."\n";}
	##
	print $OUT "$bam\t"; map { print $OUT "$counts{$_}\t" } @read_types;
	my $finished = Time::SoFar::runinterval();
	print $OUT "Time elapsed: $finished\n";
	close $OUT;
}






