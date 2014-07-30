#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/

=head1 NAME

bam2lca.pl

=head1 SYNOPSIS

Take a bam and calculate a bacterial LCA using BWA mappping against RefSeq.

=head1 DESCRIPTION


=head1 AUTHORS - Karsten Sieber

e-mail: ksieber@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

use warnings;
use strict;
use lib qw(/local/projects-t3/HLGT/scripts/lgtseek/lib/ /opt/lgtseek/lib/);      ### May need to change this depending on where the script is being run
use LGTSeek;
use Time::SoFar;
use run_cmd;
use setup_input;
use print_call;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
		'input|i=s', # Comma separated list of files
		'input_list=s',
		'Qsub=i',
		'no_gal=i',
		'hostname=s',
		'job_name=s',
		'name_sort_input=i',
		'refseq_list=s',
		'output_dir|o=s',
		'subdirs=i',
		'bin_dir=s',
		'samtools_bin=s',
		'ergatis_bin=s',
		'prinseq_bin=s',
		'donor_lineage=s',
		'host_lineage=s',
		'threads|t=i',
		'taxon_host=s',
		'taxon_dir=s',
		'taxon_idx_dir=s',
		'path_to_blastdb=s',
		'clovr=s',
		'diag=s',
		'fs=s',
		'verbose=s',
		'help|h',
		'help_full',
		'workflow_help'
		);

if($options{help}){&help;}
if(!$options{input} && !$options{input_list}){die "Error: Please give an input.bam with --input=<BAM> or --input_list=<LIST>. Try again or use --help.\n";}
## Make sure we don't end up on 
my $no_gal = undef;
my $hostname = undef;
if(defined $options{hostname}){$hostname = "$options{hostname}";$no_gal = "0";} 
else {$no_gal = defined $options{no_gal} ? "$options{no_gal}" : "1";}

## Setup Default paths for references and bins:
my $lgtseek = LGTSeek->new2(\%options);
my $inputs = setup_input(\%options);

if($options{input_list}){
	my ($fn,$dir,$suff) = fileparse(@$inputs[0],@{$lgtseek->{bam_suffix_list}});
	if($fn=~/(\w+)\_\d+$/){$fn=$1;}
	$lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}/$fn");
	$options{output_dir} = "$options{output_dir}/$fn";
}

my $i=0;
foreach my $input (@$inputs){
    $i++;	
    my ($name,$path,$suf)=fileparse($input,@{$lgtseek->{bam_suffix_list}});
	chomp $name;
	## Setup output directory
	if(!$lgtseek->{output_dir}){$lgtseek->{output_dir}="$path/lgtseq/";}
    $lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}"); 
	if($lgtseek->{subdirs}==1){
		$lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}"); 
		$lgtseek->{output_dir} = "$options{output_dir}/"."$name/"; 
		$lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}");
	}

	## Qsub the job instead of running it
	if($lgtseek->{Qsub}==1){
		my $cmd = "$^X $0";					## Start building the cmd to qsub
		if($options{input_list}){
			$options{input} = $input;										## If we are in the orignal call, we need to make sure to qsub a single input
		}
		foreach my $key (keys %options){
			next if($options{input_list} && $key=~/input_list/);			## If we are in the orignal call, we don't want to qsub more lists
			next if ($key=~/Qsub/ && !$options{input_list});			    ## If we are in the orignal call with input_list, we probably want to qsub each input
			if($options{$key}){$cmd = $cmd." --$key=$options{$key}"};		## Build the command for all other options passed in @ original call
		}
		Qsub2({
			cmd => $cmd,
			threads => "$lgtseek->{threads}",
			wd => "$lgtseek->{output_dir}",
            name => "bam2lca.$i",
            no_gal => "$no_gal",
			project => "$lgtseek->{project}",
			});
		next;
	}

	print_call(\%options);
	print STDERR "\n++++++++++++++++++++++++++++\n";
	print STDERR "+++++++  bam2lca.pl  +++++++\n";
	print STDERR "++++++++++++++++++++++++++++\n\n";

	$lgtseek->runBWA({
			input_bam => "$input",
			output_dir => "$lgtseek->{output_dir}/$name\_bwa-lca/",
			out_file => "$lgtseek->{output_dir}\/$name\_bwa-lca.txt",
			reference_list => "$lgtseek->{refseq_list}",
			threads => $lgtseek->{threads},
			run_lca => 1,
			overwrite => $lgtseek->{overwrite},
			cleanup_sai => 1,
			});

	## Cleanup Intermediate Files
    if($?==0){
        print STDERR "RM: $lgtseek->{output_dir}/$name\_bwa-lca/\n";
        $lgtseek->_run_cmd("rm -rf $lgtseek->{output_dir}/$name\_bwa-lca/");
    }
}

sub help {die "Help: This script will take bam and find the bwa-lca.
	     _____
	____/Input\\_________________________________________________________________________________
		--input=			<Input BAM>
		--input_list=			<List of BAMS> 1 per line.
	     ______
	____/Output\\________________________________________________________________________________
		--output_dir=			Directory for all output. Will be created if it doesn't exist. 
	  	  --subdirs=			<0|1> [0] 1= Make a sub-directory in output_dir based on input name
	     __________________
	____/Submit to SGE-grid\\_____________________________________________________________________
		--Qsub=				<0|1> [0] 1= qsub the job to the grid.
		  --threads=			# of threads to use for bwa sampe.
		  --job_name=			[(w+)[1..10]$] Must start with a letter character 
	_____________________________________________________________________________________________\n";
}
