#!/usr/bin/perl

=head1 NAME

lgt_search.pl

=head1 SYNOPSIS

Search an hg19.bam against human and bacteria for LGT.

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
		'input_bam=s', # Comma separated list of files
		'decrypt=s',
		'url=s',
		'split_bac_list=s',
		'hg19_ref=s',
        'refseq_list=s',
		'output_dir=s',
		'bin_dir=s',
		'samtools_bin=s',
		'ergatis_bin=s',
		'prinseq_bin=s',
		'donor_lineage=s',
		'host_lineage=s',
		'threads=s',
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
		--input_bam=				<BAM>
		--decrypt= 				[0] (0|1)
		--url=
		--split_bac_list=
		--hg19_ref=
        	--refseq_list=
		--output_dir=
		--bin_dir=
		--threads=				[1] # of CPU's to use for hyperthreading BWA. 
		--taxon_host=
		--taxon_dir=
		--taxon_idx_dir=
		--path_to_blastdb=
		--clovr=				<0|1> [0] 1=Use clovr defaults for file paths 
		--diag=					<0|1> [0] 1=Use diag node defaults for file paths 
		--fs=					<0|1> [0] 1=Use filesystem defaults for file paths 
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

my ($name,$path,$suf)=fileparse($options{input_bam},".bam");
chomp $name;

my $input_bam;
if($options{decrypt}==1){
	$input_bam = $lgtseek->decrypt({
			input => $options{input_bam},
			url => $options{url},
			output_dir => $options{output_dir}
			})
} else {
	$input_bam = $options{input_bam};
}


# Align to the donors.
print STDERR "=====RUNBWA-DONOR=====\n";
my $donor_bams = $lgtseek->runBWA({
		input_bam => $input_bam,
		output_bam => 1,
		threads => $threads,
		output_dir => "$options{output_dir}/donor_alignments/",
		reference_list => $options{split_bac_list},
		overwrite => 0,   
		cleanup_sai => 1,
});


# Align to the hosts.
print STDERR "=====RUNBWA-HOST=====\n";
my $host_bams = $lgtseek->runBWA({
		input_bam => $input_bam,
		output_bam => 1,
		threads => $threads,
		output_dir => "$options{output_dir}/host_alignments/",
		reference => $options{hg19_ref},
		overwrite => 0, 
		cleanup_sai => 1,   
});


# Postprocess the results
print STDERR "=====POSTPROCESS=====\n";
my $pp_data = $lgtseek->bwaPostProcess({
		donor_bams => $donor_bams,
		host_bams => $host_bams,
		output_prefix => $name,
		overwrite => 0,   
});

# Clean up output we don't need anymore
#print STDERR "Removing the raw donor/host mappings\n";
#print STDERR `rm -rf $options{output_dir}/host_alignments/`;
#print STDERR `rm -rf $options{output_dir}/donor_alignments/`;


# Create file with number of counts
my @header = ('run_id');
my @vals = ($name);
open OUT, ">$options{output_dir}/$name\_post_processing.tab" or die;
map {
	push(@header,$_);
	push(@vals,$pp_data->{counts}->{$_});
} ('total','host','no_map','all_map','single_map','integration_site_host','integration_site_donor','microbiome','lgt');
&print_tab("$options{output_dir}/$name\_post_processing.tab",\@header,\@vals);

# Prinseq filter the putative lgts
print STDERR "=====PRINSEQ=====\n";
my $filtered_bam = $lgtseek->prinseqFilterBam(
		{output_dir => "$options{output_dir}/prinseq_filtering",
		bam_file => $pp_data->{files}->{lgt_donor}}
		);

# Add filtered count to counts.
push(@header,'lgt_pass_prinseq_filter');
push(@vals,$filtered_bam->{count});
&print_tab("$options{output_dir}/$name\_post_processing.tab",\@header,\@vals);

sub print_tab {
	my ($file,$header,$vals) = @_;
	open OUT, ">$file" or die "Couldn't open $file\n";
	print OUT join("\t",@$header);
	print OUT "\n";
	print OUT join("\t",@$vals);
	print OUT "\n";
}

# Make LGT Fasta for Blast validation
print STDERR `mkdir -p $options{output_dir}/blast_validation`;

my $lgt_fasta = $lgtseek->sam2Fasta({
		input => "$options{output_dir}\/$name\_lgt_donor.bam"
});

# Blast & get best hits
print STDERR "=====BESTBLAST2=====\n";
my $best_blasts = $lgtseek->bestBlast2({
		db => $options{path_to_blastdb},
		lineage1 => $options{donor_lineage},
		lineage2 => $options{host_lineage},
		fasta => $lgt_fasta,
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

# Run blast and keep raw output ?
`blastall -p blastn -e 10e-5 -T F -d $options{path_to_blastdb} -i $lgt_fasta > $options{output_dir}/blast_validation/$name\_blast.raw`;

# Calculate BWA LCA's for LGTs
print STDERR "=====BWA-LCA=====\n";
print STDERR `mkdir -p $options{output_dir}/LGT_LCA-BWA`;
$lgtseek->runBWA({
		input_bam => "$options{output_dir}\/$name\_lgt_donor.bam",
		output_dir => "$options{output_dir}\/LGT_LCA-BWA/",
		out_file => "$options{output_dir}/$name\_LCA-BWA.txt",
		reference_list => "$options{refseq_list}",
		cleanup_sai => 1,
		run_lca => 1,
		overwrite => 0,
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

