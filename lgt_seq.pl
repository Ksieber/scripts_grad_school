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
my %options;
my $results = GetOptions (\%options,
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
my $threads = $options{threads} ? $options{threads} : 1;
my $lgt_coverage = $options{lgt_coverage} ? $options{lgt_coverage} : "0";
my $decrypt = $options{decrypt} ? $options{decrypt} : "0";
my $prelim_filter = $options{prelim_filter} ? $options{prelim_filter} : "0";
if($decrypt == 1 && !$options{url}){die "Error: Must give a --url to use --decrypt.\n";}

my ($name,$path,$suf)=fileparse($options{input},('.gpg.bam','_prelim.bam','.bam'));
chomp $name;
if(!$options{output_dir}){$options{output_dir}="$path/lgtseq/";}

print STDERR "+++++++++++++++++++++++++++\n";
print STDERR "++++++++  LGT-SEQ  ++++++++\n";
print STDERR "+++++++++++++++++++++++++++\n\n";
## Setup Default paths for references and bins:
my $lgtseek = LGTSeek->new2({
	options => \%options,
});


## Get input ready for processing with decryption and/or prelim filtering
print STDERR "=========PREP-INPUT========\n";
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
print STDERR "=========INPUT-READY=======\n";

__END__
# Align to the donors.
print STDERR "========RUNBWA-DONOR========\n";
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
print STDERR "========RUNBWA-HOST========\n";
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
print STDERR "=======POSTPROCESS=======\n";
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
	my $foo = $pp_data->{counts}->{$_} ? $pp_data->{counts}->{$_} : 0;
	push(@vals,$foo);
} ('total','host','no_map','all_map','single_map','integration_site_host','integration_site_donor','microbiome','lgt');
&print_tab("$options{output_dir}/$name\_post_processing.tab",\@header,\@vals);


## Check to make sure we found LGT.
print STDERR "========LGT========\n";
if($lgtseek->empty_chk({input => $pp_data->{files}->{lgt_donor}})==1){
	print STDERR "No LGT in: $pp_data->{files}->{lgt_donor}\. Skipping LGT LCA calculation and blast validation.\n";
} else {
	# Prinseq filter the putative lgts
	print STDERR "=======LGT-PRINSEQ=======\n";
	my $filtered_bam = $lgtseek->prinseqFilterBam(
		{output_dir => "$options{output_dir}/lgt_prinseq_filtering",
		input_bam => $pp_data->{files}->{lgt_host}}
	);

	# Add filtered count to counts.
	push(@header,'lgt_pass_prinseq');
	push(@vals,$filtered_bam->{count});
	&print_tab("$options{output_dir}/$name\_post_processing.tab",\@header,\@vals);

	# Calculate BWA LCA's for LGTs
	print STDERR "======LGT-BWA-LCA======\n";
	print STDERR `mkdir -p $options{output_dir}/lgt_lca-bwa`;
	$lgtseek->runBWA({
		input_bam => "$options{output_dir}\/$name\_lgt_host_filtered.bam",
		output_dir => "$options{output_dir}\/lgt_lca-bwa/",
		out_file => "$options{output_dir}\/lgt_lca-bwa\/$name\_lgt_lca-bwa.txt",
		reference_list => "$options{refseq_list}",
		run_lca => 1,
		overwrite => 0,
		cleanup_sai => 1,
	});
}

# Calculate BWA LCA's for Microbiome Reads
print STDERR "======Microbiome======\n";
## Check to make sure we found Microbiome Reads. If no microbiome reads skip this step. 
if($lgtseek->empty_chk({input => "$options{output_dir}\/$name\_microbiome.bam"})==1){
	print STDERR "No Microbiome reads in: $options{output_dir}\/$name\_microbiome.bam. Skipping microbiome LCA calculation.\n";
} else {
	# Prinseq filter the putative lgts
	print STDERR "=====Microbiome-PRINSEQ=====\n";
	my $filtered_bam = $lgtseek->prinseqFilterBam({
		output_dir => "$options{output_dir}/microbiome_prinseq_filtering",
		input_bam => "$options{output_dir}\/$name\_microbiome.bam",
	});

	# Add filtered count to counts.
	push(@header,'microbiome_pass_prinseq');
	push(@vals,$filtered_bam->{count});
	&print_tab("$options{output_dir}/$name\_post_processing.tab",\@header,\@vals);

	print STDERR "====Microbiome-BWA-LCA====\n";
	$lgtseek->runBWA({
		input_bam => "$options{output_dir}\/$name\_microbiome_filtered.bam",
		output_dir => "$options{output_dir}\/microbiome_lca-bwa/",
		out_file => "$options{output_dir}\/microbiome_lca-bwa\/$name\_microbiome_lca-bwa.txt",
		reference_list => "$options{refseq_list}",
		run_lca => 1,
		overwrite => 0,
		cleanup_sai => 1,
	}); 
}

# Calculate coverage of LGT on human side.
if($lgt_coverage==1){
	print STDERR "==Calculating Coverage of Hg19 LGT==\n";
	$lgtseek->mpileup({
		input => "$options{output_dir}\/$name\_lgt_host_filtered.bam",
		output_dir => $options{output_dir},
		ref => $options{hg19_ref},
		cleanup => 1,
		overwrite => 0,
	});
}

# Make LGT Fasta for Blast validation
my $lgt_fasta = $lgtseek->sam2Fasta({
		input => "$options{output_dir}\/$name\_lgt_host_filtered.bam",
		output_dir => "$options{output_dir}/blast_validation/"
});

# Blast & get best hits
print STDERR "=======BESTBLAST2=======\n";
my $best_blasts = $lgtseek->bestBlast2({
		db => $options{path_to_blastdb},
		lineage1 => $options{donor_lineage},
		lineage2 => $options{host_lineage},
		fasta => $lgt_fasta,
		output_dir => "$options{output_dir}/blast_validation/"
});

# Now run lgtfinder
print STDERR "========LGTFINDER=========\n";
my $valid_lgts = $lgtseek->runLgtFinder({
		lineage1 => $options{donor_lineage},
		lineage2 => $options{host_lineage},
		input_file_list => $best_blasts->{list_file},
		output_prefix => "$name",
		output_dir => "$options{output_dir}/lgt_finder/",
});

# Run blast and keep raw output ?
`blastall -p blastn -e 10e-5 -T F -d $options{path_to_blastdb} -i $lgt_fasta > $options{output_dir}/blast_validation/$name\_blast.raw`;


print STDERR "======Completed lgt_seq.pl on $options{input}======\n";


sub print_tab {
	my ($file,$header,$vals) = @_;
	open OUT, ">$file" or die "Couldn't open $file\n";
	print OUT join("\t",@$header);
	print OUT "\n";
	print OUT join("\t",@$vals);
	print OUT "\n";
}
__END__

