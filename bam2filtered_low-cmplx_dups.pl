#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/

=head1 NAME

bam_prep_for_lgt_seq.pl

=head1 SYNOPSIS

Filter a bam for low complexity reads and duplicates.

=head1 DESCRIPTION

Filter a bam for low complexity reads and duplicates.

=head1 AUTHOR - Karsten Sieber & David R. Riley

e-mail: ksieber@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut
use warnings;
use strict;
use Scalar::Util qw(reftype);
use POSIX;
use run_cmd;
use LGTSeek;
use File::Basename;
use setup_input;
use setup_output;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
my $results = GetOptions (\%options,
	'input_list=s',
	'input=s',
	'name_sort_input=i',
	'sort_mem=i',
	'Qsub=i',
	'excl=i',
	'sub_mem=i',
	'threads=i',
	'projects=s',
	'output_dir=s',
	'subdirs=i',
	'verbose=i',
	'help|h',
	);

if($options{help}){
	die "Help Basic Info: This script will remove M_M reads and split output bams into smaller bams. 
		--input=        	<BAM>
		--input_list=		<List of bams> 1 Bam / line.
		--name_sort_input= 	<0|1> [0] 1= Resort the input bam by read names.  
		--tmp_dir=   		Directory for intermediate files. Default = input_dir/prinseq/
		--output_dir=
		--Qsub=
		  --excl=
		  --sub_mem=
		  --threads=\n";
}

my $lgtseek = LGTSeek->new2(\%options);
$options{name_sort_input} = defined $options{name_sort_input} ? $options{name_sort_input} : "0";
my $input = setup_input(\%options);	## $input is a array ref. 

my $i=0;
foreach my $input (@$input){
	$i++;
    my ($fn,$path,$suf)=fileparse($input,@{$lgtseek->{bam_suffix_list}});
	$options{output_dir} = defined $options{output_dir} ? $options{output_dir} : $path;
	my $output_dir = defined $options{subdirs} ? "$options{output_dir}/$fn/" : $options{output_dir};
	$lgtseek->_run_cmd("mkdir -p $output_dir");
	if($lgtseek->{Qsub}==1){
		## If we are in the orignal call, change input from list to a single file
		if($options{input_list}){$options{input} = $input;}
		## Check $sub_mem is enough for sorting
		my $original_sub_mem; 
		my $original_sort_mem;
		if($options{sub_mem} =~ /^(\d+)[K|M|G]$/){$original_sub_mem = $1;}
		if($options{sort_mem} =~ /^(\d+)[K|M|G]$/){$original_sort_mem = $1;}
		if($original_sub_mem < ($original_sort_mem * $options{threads})){
			$options{sub_mem} = (ceil(($original_sort_mem * $options{threads})* 1.1) )+1 . "G";
		}
		## Build qsub command
		my $cmd = "/home/ksieber/scripts/bam2filtered_low-cmplx_dups.pl"; 
		foreach my $key (keys %options){
			next if ($key=~/Qsub/);
			next if($options{input_list} && $key=~/input_list/);			  ## If we are in the orignal call, we don't want to qsub more lists
			if(defined $options{$key}){$cmd = $cmd." --$key=$options{$key}"};
		}
	    ## submit command to grid
	    Qsub2({
	    	cmd => "$cmd",
	    	wd => "$output_dir",
	    	name => "bam2filter",
	    	mem => "$lgtseek->{sub_mem}",
	    	threads => "$lgtseek->{threads}",
	    	project => "$lgtseek->{project}",
	    	excl => "$options{excl}",
	    	});
		## Skip to next input for qsub
		next;
	}

	print STDERR "======== bam2filtered: Start ========\n";
	## name sort input
	if($options{name_sort_input}==1){
		print STDERR "======== Name Sorting ========\n";
		my ($fn,$path,$suf)=fileparse($input,@{$lgtseek->{bam_suffix_list}});
		run_cmd("samtools sort -n -@ $lgtseek->{threads} -m $lgtseek->{sort_mem} $input $output_dir/$fn\.name-sort");
		$input = "$output_dir/$fn\.name-sort.bam";
	}

	## PrinSeq Filter
	my $filtered_bams = $lgtseek->prinseqFilterBam({
		input_bam => $input,
		output_dir => $output_dir,
		overwrite => $lgtseek->{overwrite},
		});
}






										


