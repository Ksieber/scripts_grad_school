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
	'sort_mem=s',
	'split_bam=s',
	'encrypt=s',
	'key=s',
	'prelim_filter=s',
	'seqs_per_file=s',
	'keep_softclip=s',
	'Qsub=s',
	'sub_mem=s',
	'projects=s',
	'output_dir=s',
	'subdirs=s',
	'output_list',
	'overwrite=s',
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
		----------------------------------------------------------------------------------------
		--input=		<BAM>
		--input_list=   	<LIST of BAMS> 1 bam per line.
		----------------------------------------------------------------------------------------
		--prelim_filter= 	<0|1> [1] 1= Filter out M_M reads, keeping MU,UU,and SC. 
		  --keep_softclip=	<0|1> [1] 1= Keep soft clipped reads >=24 bp (Reads potentially on LGT) 
		--name_sort_input= 	<0|1> [1] 1= Resort the input bam by read names.  
		  --sort_mem=		<Bytes> [7500000000] Mem free to sort reads with. More memory limits intermidiate files.  
		--split_bam=		<0|1> [1] 1= Split bam(s)
		  --seqs_per_file=    	<lines per bam> [50,000,000]
		--encrypt=       	<0|1> [0] 1= encrypt ** untested **
		  --key=          	GPG key to use for encryption.
		----------------------------------------------------------------------------------------
		--Qsub=          	<0|1> [0] 1= Qsub this script for each input. 
		  --project=		<project> [jdhotopp-lab] SGE group project to submit command under.
		  --sub_mem=		<mem min> [40G] Min memory available at submission. ** Make sure enough memory for --name_sort_input
		----------------------------------------------------------------------------------------
		--output_dir=   	Directory for output. 
		  --subdirs=       	<0|1> [0] 1= Make a directory under the output_dir to place the output. 
		  --output_list=  	<0|1> [1] 1= Make a list of the output created.
		----------------------------------------------------------------------------------------
		--overwrite=		<0|1> [0] 1= Overwrite previous files.
		----------------------------------------------------------------------------------------
		--bin_dir=     		[/local/projects-t3/HLGT/scripts/lgtseek/bin/] Directory where LGTSeek.pm is stored.
		--fs=			<0|1> [1] 1= IGS filesystem.
		--clovr=		<0|1> [0] 1= clovr box
		--diag=			<0|1> [0] 1= Diag node
		----------------------------------------------------------------------------------------
		--help 			Basic Help Info
		--help_full		Full Help Info
		----------------------------------------------------------------------------------------\n";
}

if(!$options{input} && !$options{input_list}){die "Error: Must give input with --input=<BAM> or --input_list=<LIST of bams>\n";}

## Set default values
$options{name_sort_input} = $options{name_sort_input} ? $options{name_sort_input} : 1;
$options{prelim_filter} = $options{prelim_filter} ? $options{prelim_filter} : 1;
$options{keep_softclip} = $options{keep_softclip} ? $options{keep_softclip} : 1;
$options{split_bam} = $options{split_bam} ? $options{split_bam} : 1;
$options{overwrite} = $options{overwrite} ? $options{overwrite} : 0;
$options{sub_mem} = $options{sub_mem} ? $options{sub_mem} : "10G";

my $lgtseek = LGTSeek->new2({
	options => \%options,											## Pass %options to lgtseek.
	});

my $input = setup_input(\%options);	## $input is a array ref. 

foreach my $input (@$input){
	my ($fn,$path,$suf)=fileparse($input,('_resorted.bam','.bam'));
	my $output_dir = $options{subdirs} ? "$options{output_dir}/$fn/" : $options{output_dir};
	$lgtseek->_run_cmd("mkdir -p $output_dir");
	
	## Qsub this script foreach input and any of the options passed
	if($lgtseek->{Qsub}==1){
		my $cmd = "/home/ksieber/scripts/lgt_prep_bam_lgt_seq.pl"; 
		foreach my $key (keys %options){
			next if ($key=~/Qsub/);
			if($options{$key}){$cmd = $cmd." --$key=$options{$key}"};
		}
        $fn =~ /(\w{1,10})$/;
        my $job_name = $1;
		Qsub2({
			cmd => $cmd,
			wd => $output_dir,
			name => $job_name,
			mem => $options{sub_mem},
            project => $lgtseek->{project},
			});
		next;
	}
	
	## PrelimFilter to remove M_M reads. 
	my $bams;
	if($lgtseek->{prelim_filter}==1 || $lgtseek->{split_bam}==1 || $lgtseek->{name_sort_input}==1){
		$bams = $lgtseek->prelim_filter({
			input_bam => $input,
			output_dir => $output_dir,
			name_sort_input => $lgtseek->{name_sort_input},				## Default = 1
			sort_mem => $lgtseek->{sort_mem},							## Default = 6G lgtseek default. lgt_prep overides to 40G. 
			split_bam => $lgtseek->{split_bam},							## Default = 1
			seqs_per_file => $lgtseek->{seqs_per_file},					## Default = 50M
			keep_softclip => $lgtseek->{keep_softclip},					## Default = 1
			overwrite => $lgtseek->{overwrite},							## Default = 0
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
		open(my $olistfh, ">$output_dir/output.list") or die "Unable to open: $output_dir/output.list because: $!\n";
 		if($lgtseek->{encrypt}==1){foreach my $out (@encrypted){print $olistfh "$out\n";}} 
		elsif ($lgtseek->{split_bam}==1 || $lgtseek->{prelim_filter}==1){foreach my $out2 (@$bams){print $olistfh "$out2\n";}}
		#else {print $olistfh "$in_bam\n";}
	}
}


__END__

