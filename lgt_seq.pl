#!/usr/bin/perl

=head1 NAME

lgt_search.pl

=head1 SYNOPSIS

Search an hg19.bam against human and bacteria for LGT.

=head1 DESCRIPTION


=head1 AUTHORS - Karsten Sieber & David Riley

e-mail: ksieber@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut
use warnings;
use strict;
use lib qw(/local/projects-t3/HLGT/scripts/lgtseek/lib/ /opt/lgtseek/lib/);      ### May need to change this depending on where the script is being run
use LGTSeek;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use setup_default_paths;
our %options;
our $results = GetOptions (\%options,
		'input=s', # Comma separated list of files
		'decrypt=s',
		'url=s',
		'prelim_filter=s',
		'split_bac_list=s',
		'hg19_ref=s',
        'refseq_list=s',
		'output_dir=s',
		'lgt_coverage=s',
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
		--input=				<BAM>
		--decrypt= 				<0|1> [0]
		--url=
		--prelim_filter			<0|1> [0] 1=Filter out human M_M reads from original input.
		--split_bac_list=
		--hg19_ref=
        	--refseq_list=
		--output_dir=
		--lgt_coverage=				<0|1> [0] 1= Calculate coverage of hg19 LGT. 
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

if(!$options{input}){die "Error: Please give an input.bam with --input=<FILE>. Try again or use --help.\n";}
if(!$options{output_dir}){print "It is HIGHLY recommended you STOP, restart, and use a --output_dir=<some/where/>.\n";sleep 60;}



# Take care of the inputs
## Setup Default paths for references and bins:
&setup_default_paths();

## These are only set to : <X> if --options{y} isn't used OR --clovr/diag/fs=1
my $bin_dir = $options{bin_dir} ? $options{bin_dir} : '/opt/lgtseek/bin/';    
my $ergatis_dir = $options{ergatis_bin} ? $options{ergatis_bin} :'/opt/ergatis/bin/';
my $prinseq_bin = $options{prinseq_bin} ? $options{prinseq_bin} : '/opt/prinseq/bin/';
my $samtools_bin = $options{samtools_bin} ? $options{samtools_bin} : 'samtools';
my $threads = $options{threads} ? $options{threads} : 1;
my $lgt_coverage = $options{lgt_coverage} ? $options{lgt_coverage} : "0";
my $decrypt = $options{decrypt} ? $options{decrypt} : "0";
my $prelim_filter = $options{prelim_filter} ? $options{prelim_filter} : "0";
if($decrypt == 1 && !$options{url}){die "Error: Must give a --url to use --decrypt.\n";}

my ($name,$path,$suf)=fileparse($options{input},('.gpg.bam','_prelim.bam','.bam'));
chomp $name;
if(!$options{output_dir}){$options{output_dir}="$path/lgtseq/";}

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

## Get input ready for processing with decryption and/or prelim filtering
print STDERR "=====PREP-INPUT=====\n";
my $input_bam;
if($decrypt==1){
	$input_bam = $lgtseek->decrypt({
			input => $options{input},
			url => $options{url},
			output_dir => $options{output_dir}
	})
} elsif($prelim_filter==1){
	my $unfiltered_bam;
	if ($decrypt==1) {$unfiltered_bam = $input_bam;} 
	else {$unfiltered_bam = $options{input};}
	$input_bam = $lgtseek->prelim_filter({
			input_bam => $unfiltered_bam,
			output_dir => "$options{output_dir}/prelim_filter/",
			keep_softclip => 1,
			overwrite => 0,
	});
} else {$input_bam = $options{input};}



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

if($lgtseek->empty_chk({input => $pp_data->{files}->{lgt_donor}})==1){die "No LGT in: $pp_data->{files}->{lgt_donor}\nStopping lgt_seq.\n"}

# Prinseq filter the putative lgts
print STDERR "=====PRINSEQ=====\n";
my $filtered_bam = $lgtseek->prinseqFilterBam(
		{output_dir => "$options{output_dir}/prinseq_filtering",
		input_bam => $pp_data->{files}->{lgt_donor}}
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
my $lgt_fasta = $lgtseek->sam2Fasta({
		input => "$options{output_dir}\/$name\_lgt_donor.bam",
		output_dir => "$options{output_dir}/blast_validation/"
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
print STDERR "=====LGT-BWA-LCA=====\n";
print STDERR `mkdir -p $options{output_dir}/lgt_lca-bwa`;
$lgtseek->runBWA({
		input_bam => "$options{output_dir}\/$name\_lgt_donor.bam",
		output_dir => "$options{output_dir}\/lgt_lca-bwa/",
		out_file => "$options{output_dir}\/lgt_lca-bwa\/$name\_lgt_lca-bwa.txt",
		reference_list => "$options{refseq_list}",
		run_lca => 1,
		overwrite => 0,
		cleanup_sai => 1,
}); 

# Calculate BWA LCA's for Microbiome Reads
print STDERR "=====Microbiome-BWA-LCA=====\n";
print STDERR `mkdir -p $options{output_dir}/microbiome_lca-bwa`;
$lgtseek->runBWA({
		input_bam => "$options{output_dir}\/$name\_microbiome.bam",
		output_dir => "$options{output_dir}\/microbiome_lca-bwa/",
		out_file => "$options{output_dir}\/microbiome_lca-bwa\/$name\_lgt_lca-bwa.txt",
		reference_list => "$options{refseq_list}",
		run_lca => 1,
		overwrite => 0,
		cleanup_sai => 1,
}); 

if($lgt_coverage==1){
	print STDERR "=====Calculating Coverage of Hg19 LGT======\n";
	$lgtseek->mpileup({
		input => "$options{output_dir}\/$name\_lgt_donor.bam",
		output_dir => $options{output_dir},
		ref => $options{hg19_ref},
		cleanup => 1,
		overwrite => 0,
	});
}

print STDERR "Completed lgt_seq.pl on $options{input}\n";

__END__

