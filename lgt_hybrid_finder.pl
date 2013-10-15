#!/usr/bin/perl

=head1 NAME

lgt_hybrid_finder.pl

=head1 SYNOPSIS

Search for reads ON the integration point. ie Half Donor Half Host.

=head1 DESCRIPTION


=head1 AUTHOR - Karsten B. Sieber

e-mail: ksieber@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

use strict;
use lib '/local/projects-t3/HLGT/scripts/lgtseek/lib/';      ### May need to change this depending on where the script is being run
use LGTSeek;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
		'input=s', # Comma separated list of files
		'output_dir=s',
		'bin_dir=s',
		'samtools_bin=s',
		'ergatis_bin=s',
		'prinseq_bin=s',
		'donor_lineage=s',
		'host_lineage=s',
		'taxon_host=s',
		'taxon_dir=s',
		'taxon_idx_dir=s',
		'path_to_blastdb=s',
		'clovr=s',
		'diag=s',
		'fs=s',
		'help|h'
		);

if($options{help}){die "Help: This script will takes a bam and identifies bacterial human LGT.
		--input=				<.bam or .fa>
		--output_dir=				</dir/for/output/> [cwd]
		--taxon_host=
		--taxon_dir=
		--taxon_idx_dir=
		--path_to_blastdb=
		--clovr=				<0|1> [0] 1=Use clovr defaults for file paths 
		--diag=					<0|1> [0] 1=Use diag node defaults for file paths 
		--fs=					<0|1> [1] 1=Use filesystem defaults for file paths 
		--bin_dir=
		--samtools_bin=
		--ergatis_bin=
		--help\n";
}

if(!$options{input_bam}){die "Error: Please give an input.bam with --input_bam=<FILE>. Try again or use --help.\n";}
if($options{decrypt} == 1 && !$options{url}){die "Error: Must give a --url to use --decrypt.\n";}

# Take care of the inputs
## Setup Default paths for references and bins:
setup_defaults();
my $bin_dir = $options{bin_dir} ? $options{bin_dir} : '/opt/lgtseek/bin/';    
my $ergatis_dir = $options{ergatis_bin} ? $options{ergatis_bin} :'/opt/ergatis/bin/';
my $prinseq_bin = $options{prinseq_bin} ? $options{prinseq_bin} : '/opt/prinseq/bin/';
my $samtools_bin = $options{samtools_bin} ? $options{samtools_bin} : 'samtools';
my $threads = $options{threads} ? $options{threads} : 1;
my @in_suffix_list=(".bam",".fa");

# Create an lgtseek object
my $lgtseek = LGTSeek->new({
		bin_dir => $bin_dir,
		output_dir => $options{output_dir},
		ergatis_bin => $ergatis_dir,
		prinseq_bin => $prinseq_bin,
		paired_end => 1,
		taxon_host => $options{taxon_host},
		taxon_dir => $options{taxon_dir},
		taxon_idx_dir => $options{taxon_idx_dir}
});

my ($name,$path,$suf)=fileparse($options{input},@in_suffix_list);
chomp $name;
my $input;
if($suf=~/bam/){
	$input = $lgtseek->sam2Fasta({input => "$options{input}"});
} elsif ($suf=~/fa/){
	$input=$options{input};
} else {die "Error: Could not determine input file type. Please use a .bam or .fa\n";}

# Blast & get best hits
print STDERR "=====BESTBLAST2=====\n";
my $best_blasts = $lgtseek->bestBlast2({
		db => $options{path_to_blastdb},
		lineage1 => $options{donor_lineage},
		lineage2 => $options{host_lineage},
		fasta => $input,
		output_dir => "$options{output_dir}/blast_validation/"
});

# Now run lgtfinder
print STDERR "=====LGTFINDER=====\n";
my $valid_lgts = $lgtseek->runLgtFinder({
		lineage1 => $options{donor_lineage},
		lineage2 => $options{host_lineage},
		input_file_list => $best_blasts->{list_file},
		output_prefix => "$name",
		output_dir => "$options{output_dir}/lgt_finder/",
});


sub setup_defaults {
## Default file Path prefixes
	my $diag = 
	{
		bin_dir => "/opt/lgtseek/bin/",
		ergatis_bin => "/opt/ergatis/bin/",
		prinseq_bin => "/opt/prinseq/bin/",
		samtools_bin => "samtools",
		split_bac_list => "/mnt/staging/data/lgt_seq/mnt/references/split_refseq_bacteria/split_bacteria_ref.list",
		hg19_ref => "/mnt/staging/data/lgt_seq/mnt/references/hg19/dna/hg19.fa",
		refseq_list => "/mnt/staging/data/lgt_seq/mnt/references/refseq_bacteria_BWA_INDEXED_20110831/refseq.list",
		taxon_host => "cloud-128-152.diagcomputing.org:10001",
		taxon_dir => "/mnt/staging/data/lgt_seq/mnt/references/taxonomy",
		taxon_idx_dir => "/mnt/staging/data/lgt_seq/mnt/references/taxonomy",
		path_to_blastdb => "/mnt/scratch/ksieber/ref/refseq_bacteria_merged.fna",  ## FIX this
		donor_lineage => "Bacteria",
		host_lineage => "Eukaryota"
	};
	
	my $clovr = 
	{
		bin_dir => "/opt/lgtseek/bin/",
		ergatis_bin => "/opt/ergatis/bin/",
		prinseq_bin => "/opt/prinseq/bin/",
		samtools_bin => "samtools",
		split_bac_list => "/mnt/references/split_refseq_bacteria/split_bacteria_ref.list",
		hg19_ref => "/mnt/references/hg19/dna/hg19.fa",
		refseq_list => "/mnt/references/refseq_bacteria_BWA_INDEXED_20110831/refseq.list",
		taxon_host => "cloud-128-152.diagcomputing.org:10001",						## FIX this
		taxon_dir => "/mnt/references/taxonomy/20120720_135302/",
		taxon_idx_dir => "/mnt/references/taxonomy/20120720_135302/",
		path_to_blastdb => "/mnt/scratch/ksieber/ref/refseq_bacteria_merged.fna",  ## FIX this
		donor_lineage => "Bacteria",
		host_lineage => "Eukaryota"
	};

	my $fs = 
	{
		bin_dir => "/local/projects-t3/HLGT/scripts/lgtseek/bin/",
		ergatis_bin => "/local/projects/ergatis/package-driley/bin/",
		prinseq_bin => "/home/ksieber/lib/prinseq-lite-0.18.1/",
		samtools_bin => "samtools",
		split_bac_list => "/local/projects-t3/HLGT/references/split_bacteria/all_bacteria.list",
		hg19_ref => "/local/projects-t3/HLGT/references/hg19/hg19.fa",
		refseq_list => "/local/projects-t3/HLGT/references/refseq_bacteria_BWA_INDEXED_20110831/refseq.list",
		taxon_host => "mongotest1-lx.igs.umaryland.edu:10001",
		taxon_dir => "/local/db/repository/ncbi/blast/20120414_001321/taxonomy/",
		taxon_idx_dir => "/local/projects-t3/HLGT/idx_dir/20120414",
		path_to_blastdb => "/local/db/ncbi/blast/db/nt",  ## FIX this
		donor_lineage => "Bacteria",
		host_lineage => "Eukaryota"
	};
	
	if($options{clovr}==1){foreach my $keys (keys %$clovr){$options{$keys}=$clovr->{$keys};}}
	if($options{diag}==1){foreach my $keys (keys %$diag){$options{$keys}=$diag->{$keys};}}			
	if($options{fs}==1){foreach my $keys (keys %$fs){$options{$keys}=$fs->{$keys};}}
}
