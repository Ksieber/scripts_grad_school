#!/usr/bin/perl

=head1 NAME

bam_prep_for_lgt_seq.pl

=head1 SYNOPSIS

Split, encrypt, and/or prelim filter a bam for lgt_seq.pl 

=head1 DESCRIPTION

Split, encrypt, and/or prelim filter a bam for lgt_seq.pl 

=head1 AUTHOR - Karsten Sieber & David R. Riley

e-mail: ksieber@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut
#use warnings;
use strict;
use Scalar::Util qw(reftype);
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
	'name_sort_input=s',
	'split_bam=s',
	'encrypt=s',
	'key=s',
	'prelim_filter=s',
	'seqs_per_file=s',
	'keep_softclip=s',
	'Qsub=s',
	'output_dir=s',
	'subdirs=s',
	'output_list',
	'samtools_bin=s',
	'ergatis_dir=s',
	'output_list=s',
	'bin_dir=s',
	'fs=s',
	'clovr=s',
	'diag',
	'help|h',
	'help_full'
	);

if($options{help}){
	die "Help Basic Info: This script will remove M_M reads and split output bams into smaller bams. 
		--input=        	<BAM>
		--name_sort_input= 	<0|1> [0] 1= Resort the input bam by read names.  
		--output_dir=   	Directory for output. 
		--help 			Basic Help Info
		--help_full		Full Help Info\n";
}

if($options{help_full}){
	die "Help Full Info: This script will remove M_M reads, keeping M_U, U_U, M_M with Softclip. The output bams are split into smaller chunks. 
		--input=			<BAM>
		--input_list=   	<LIST of BAMS> 1 bam per line.
		--name_sort_input= 	<0|1> [1] 1= Resort the input bam by read names.  
		--prelim_filter= 	<0|1> [1] 1= Filter out M_M reads, keeping MU,UU,and SC. 
		--keep_softclip=	<0|1> [1] 1= Keep soft clipped reads >=24 bp (Reads potentially on LGT) 
		--split_bam=		<0|1> [1] 1= Split bam(s)
		--seqs_per_file=    	<lines per bam> [50,000,000]
		--encrypt=       	<0|1> [0] 1= encrypt ** untested **
		--key=          	GPG key to use for encryption.
		--Qsub=          	<0|1> [0] 1= Qsub this script for each input. 
		--output_dir=   	Directory for output. 
		--subdirs=       	<0|1> [0] 1= Make a directory within the output_dir to place the output. 
		--output_list=  	<0|1> [0] 1= Make a list of the output created.
		--bin_dir=      	Directory where LGTSeek.pm is stored.
		--fs=			<0|1> [1] 1= IGS filesystem.
		--clovr=		<0|1> [0] 1= clovr box
		--diag=			<0|1> [0] 1= Diag node
		--help 			Basic Help Info
		--help_full		Full Help Info\n";
}

if(!$options{input} && !$options{input_list}){die "Error: Must give input with --input=<BAM> or --input_list=<LIST of bams>\n";}
$options{name_sort_input} = $options{name_sort_input} ? $options{name_sort_input} : 1;
$options{prelim_filter} = $options{prelim_filter} ? $options{prelim_filter} : 1;
$options{keep_softclip} = $options{keep_softclip} ? $options{keep_softclip} : 1;
$options{split_bam} = $options{split_bam} ? $options{split_bam} : 1;

my $lgtseek = LGTSeek->new2({
	options => \%options,
	});

my $input = setup_input(\%options);	## $input is a array ref. 

foreach my $input (@$input){
	my ($fn,$path,$suf)=fileparse($input,('_resorted.bam','.bam'));
	## Qsub this script foreach input and any of the options passed
	if($lgtseek->{Qsub}==1){
		my $output_dir = $options{subdirs} ? "$options{output_dir}/$fn/" : $options{output_dir};
		$lgtseek->_run_cmd("mkdir -p $output_dir");
		$lgtseek->_run_cmd('cd $output_dir');
		my $cmd = "/home/ksieber/scripts/bam_prep_for_lgt_seq.pl";
		foreach my $key (keys %options){
			next if ($key=~/Qsub/);
			if($options{$key}){$cmd = $cmd." --$key=$options{$key}"};
		}
		Qsub2({
			cmd => $cmd,
			threads => "$options{threads}",
			wd => $output_dir,
			});
		last;
	}
	my $out_dir;
	if($lgtseek->{subdirs}==1){
		$out_dir = "$lgtseek->{output_dir}/$fn/";
	} else {
		$out_dir = "$lgtseek->{output_dir}";
	}

	## PrelimFilter to remove M_M reads. 
	my $bams;
	if($lgtseek->{prelim_filter}==1 || $lgtseek->{split_bam}==1 || $lgtseek->{name_sort_input}==1){
		$bams = $lgtseek->prelim_filter({
			input_bam => $input,
			output_dir => $out_dir,
			});
	}

	## Encrypt the file(s) 
	my @encrypted;
	if($lgtseek->{encrypt}==1){
		foreach my $files (@$bams){
			my $cmd = "gpg -o $files\.gpg --default-key $options{key} -r $options{key} -e $files";
			if ($lgtseek->{Qsub}==1){Qsub($cmd);} 
			else {$lgtseek->_run_cmd($cmd);}
			push(@encrypted,"$files\.gpg");
		}
	}

	## Print out the output list
	if($lgtseek->{output_list}==1){
		# Open a list file to write the output bams to
		my $olistfh;
		if($lgtseek->{output_list}==1) {
			open($olistfh, ">$out_dir/output.list") or die "Unable to open: $out_dir/output.list because: $1\n";
		} else {
			$olistfh = *STDERR;
		}
 		if($lgtseek->{encrypt}==1){foreach my $out (@encrypted){print $olistfh "$out\n";}} 
		elsif ($lgtseek->{split_bam}==1 || $lgtseek->{prelim_filter}==1){foreach my $out2 (@$bams){print $olistfh "$out2\n";}}
		#else {print $olistfh "$in_bam\n";}
	}
}


__END__

