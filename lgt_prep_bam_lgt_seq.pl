#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/

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
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
my $results = GetOptions (\%options,
	'input_list=s',
	'input|i=s',
	'name_sort_input=s',
	'sort_mem=s',
	'split_bam=i',
	'encrypt=i',
	'key=s',
	'prelim_filter=i',
	'seqs_per_file=i',
	'keep_softclip=i',
	'Qsub=i',
	'excl=i',
	'sub_mem=s',
	'threads|t=i',
	'projects=s',
	'output_dir|o=s',
	'subdirs=i',
	'overwrite=s',
	'samtools_bin=s',
	'ergatis_dir=s',
	'output_list=s',
	'bin_dir=s',
	'fs=s',
	'clovr=s',
	'diag',
	'verbose=i',
	'print_hostname|ph=i',
	'config_file=s',
	'help|h',
	'help_full'
	);
use print_call;
print_hostname(\%options);									## This is useful for trouble shooting grid nodes that might be missing modules for LGTSeek etc. 
use Scalar::Util qw(reftype);
use POSIX;
use run_cmd;
use LGTSeek;
use File::Basename;
use setup_input;
use setup_output;

if($options{help}){
	die "Help Basic Info: This script will remove M_M reads and split output bams into smaller bams. 
		--input=        	<BAM>
		--name_sort_input= 	<0|1> [0] 1= Resort the input bam by read names.  
		--output_dir=   	Directory for output. 
		--help_full		Full Help Info\n";
}

if($options{help_full}){
	die "Help Full Info: This script will remove M_M reads, keeping M_U, U_U, M_M with Softclip. The output bams are split into smaller chunks. 
		----------------------------------------------------------------------------------------
		--input=		<BAM>
		--input_list=   	<LIST of BAMS> 1 bam per line.
		----------------------------------------------------------------------------------------
		--paired_end=		<0|1> [1] 1= Paired End Sequencing reads.
		--prelim_filter= 	<0|1> [1] 1= Filter out M_M reads, keeping MU,UU,and SC. 
		  --keep_softclip=	<0|1> [1] 1= Keep soft clipped reads >=24 bp (Reads potentially on LGT) 
		--name_sort_input= 	<0|1> [1] 1= Resort the input bam by read names.  
		  --sort_mem=		[5G] Mem per thread to sort with. 
		  --threads=		[1] # of threads 
		--split_bam=		<0|1> [1] 1= Split bam(s)
		  --seqs_per_file=    	<lines per bam> [50,000,000]
		--encrypt=       	<0|1> [0] 1= encrypt ** untested **
		  --key=          	GPG key to use for encryption. ** untested **
		----------------------------------------------------------------------------------------
		--Qsub=          	<0|1> [0] 1= Qsub this script for each input. 
		  --project=		<project> [jdhotopp-lab] SGE group project to submit command under.
		  --sub_mem=		[7G] --sub_mem MUST >= (--threads * --sort_mem)
		  --hostname=		Designate a specific host to run on
		  --no_gal=		<0|1> [1] 1= Specify not to run on galactus b/c lgtseek fails w/ LGTSeek.pm
		  --excl=		<0|1> [0] 1= Run exclusively on a node.
		----------------------------------------------------------------------------------------
		--output_dir=   	Directory for output. 
		  --subdirs=       	<0|1> [0] 1= Make a directory under the output_dir to place the output. 
		  --output_list=  	<0|1> [1] 1= Make a list of the output created.
		----------------------------------------------------------------------------------------
		--overwrite=		<0|1> [0] 1= Overwrite previous files.
		----------------------------------------------------------------------------------------
		--bin_dir=     		[/local/projects-t3/HLGT/scripts/lgtseek/bin/] Directory where LGTSeek.pm is stored.
		----------------------------------------------------------------------------------------
		--verbose		<0|1> [1] 1= Verbose reporting of progress. 0 =Turns off reports. 
		--help 			Basic Help Info
		--help_full		Full Help Info
		--conf_file=				[~/.lgtseek.conf]
		----------------------------------------------------------------------------------------\n";
}

if(!$options{input} && !$options{input_list}){die "Error: Must give input with --input=<BAM> or --input_list=<LIST of bams>\n";}

## Set default values
$options{name_sort_input} = defined $options{name_sort_input} ? "$options{name_sort_input}" : "1";
$options{prelim_filter}   = defined $options{prelim_filter}   ? "$options{prelim_filter}"   : "1";
$options{keep_softclip}   = defined $options{keep_softclip}   ? "$options{keep_softclip}"   : "1";
$options{split_bam}       = defined $options{split_bam}       ? "$options{split_bam}"       : "1";
$options{overwrite}       = defined $options{overwrite}       ? "$options{overwrite}"       : "0";
$options{sub_mem}         = defined $options{sub_mem}         ? "$options{sub_mem}"         : "7G";
$options{threads}         = defined $options{threads}         ? "$options{threads}"         : "1";
$options{excl} 			  = defined $options{excl} 			  ? "$options{excl}" 			: "0";
my $no_gal = undef;
my $hostname = undef;
if(defined $options{hostname}){
	$hostname = "$options{hostname}";
	$no_gal = "0";
} else {
	$no_gal = defined $options{no_gal} ? "$options{no_gal}" : "1";
}


my $lgtseek = LGTSeek->new2(\%options);
my $input = setup_input(\%options);	## $input is a array ref. 

foreach my $input (@$input){
	my ($fn,$path,$suf)=fileparse($input,('_resorted.bam','.bam'));
	my $output_dir = $options{subdirs} ? "$options{output_dir}/$fn/" : $options{output_dir};
	$lgtseek->_run_cmd("mkdir -p $output_dir");
	
	## Qsub this script foreach input and any of the options passed
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
		my $cmd = "/home/ksieber/scripts/lgt_prep_bam_lgt_seq.pl"; 
		foreach my $key (keys %options){
			next if ($key=~/Qsub/);
			next if($options{input_list} && $key=~/input_list/);			  ## If we are in the orignal call, we don't want to qsub more lists
			if(defined $options{$key}){$cmd = $cmd." --$key=$options{$key}"};
		}
        $fn =~ /(\w{1,10})$/;
        my $job_name = $1;

        ## submit command to grid
		Qsub2({
			cmd => "$cmd",
			wd => "$output_dir",
			name => "$job_name",
			mem => "$options{sub_mem}",
			threads => "$lgtseek->{threads}",
            project => "$lgtseek->{project}",
            hostname => "$hostname",
            no_gal => "$no_gal",
            excl => "$options{excl}",
			});

		## Skip to next input for qsub
		next;									
	}
	
	## PrelimFilter to remove M_M reads. 
	my $bams;
	if($lgtseek->{prelim_filter}==1 || $lgtseek->{split_bam}==1 || $lgtseek->{name_sort_input}==1){
		$bams = $lgtseek->prelim_filter({
			input_bam => $input,
			output_dir => $output_dir,
			name_sort_input => $lgtseek->{name_sort_input},				## Default = 1
			sort_mem => $lgtseek->{sort_mem},							## Default = 10G lgtseek default. lgt_prep overides to 40G. 
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
	}
	print STDERR "Completed lgt_prep_bam_lgt_seq.pl on: $input\n";
}


__END__

