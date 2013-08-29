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
use Time::SoFar;
use run_cmd;
use setup_input;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
		'input=s', # Comma separated list of files
		'input_list=s',
		'Qsub=s',
		'decrypt=s',
		'url=s',
		'prelim_filter=s',
		'name_sort_input=s',
		'keep_softclip=s',
		'split_bam=s',
		'seqs_per_file=s',
		'split_bac_list=s',
		'hg19_ref=s',
		'refseq_list=s',
		'output_dir=s',
		'subdirs=s',
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
		'help|h',
		'help_full'
		);

if($options{help}){die "Help: This script will identify bacterial human LGT.
		--input=			<Input BAM>
		--output_dir=			Directory for all output. Will be created if it doesn't exist. 
		--subdirs=			<0|1> [0] 1= Make a sub-directory in output_dir based on input name
		--help				Help Basic Info
		--help_full 			Help Full Info\n";
}

if($options{help_full}){die "Help: This script takes a bam and identifies bacterial human LGT.
		--input=			<Input BAM>
		--input_list=		<List of BAMS> 1 per line.
		--decrypt= 			<0|1> [0] 1= Decrypt the input bam with a key downloaded from --url=
		--url=				The url to download the decryption key from.
		--name_sort_input=		<0|1> [0] 1= Resort input bam by read name. 
		--prelim_filter=		<0|1> [0] 1= Filter out human M_M reads from original input.
		--keep_softclip=		<0|1> [1] 1= Keep reads that are softclipped >=24 bp. 
		--split_bam=			<0|1> [1] 1= Split bam into --seqs_per_file chunks. 
		--seqs_per_file=		[50000000] (50 Million)
		--split_bac_list=		Path to the list of split bacterial references (A2D, E2P, R2Z)
		--hg19_ref=			Path to hg19 reference
		--refseq_list=			Path to all bacterial references in refseq. 
		--output_dir=			Directory for all output. Will be created if it doesn't exist. 
		--subdirs=			<0|1> [0] 1= Make a sub-directory in output_dir based on input name
		--Qsub=				<0|1> [0] 1= qsub the job to the grid.
		--project=			[jdhotopp-lab] Grid project to use. 
		--lgt_coverage=			<0|1> [0] 1= Calculate coverage of hg19 LGT. 
		--bin_dir=			[/local/projects-t3/HLGT/scripts/lgtseek/bin/]
		--threads=			[1] # of CPU's to use for multithreading BWA sampe 
		--taxon_host=			[mongotest1-lx.igs.umaryland.edu:10001]
		--taxon_dir=			[/local/db/repository/ncbi/blast/20120414_001321/taxonomy/]
		--taxon_idx_dir=		[/local/projects-t3/HLGT/idx_dir/20120414]
		--path_to_blastdb=		[/local/db/repository/ncbi/blast/20120414_001321/nt/nt]
		--clovr=			<0|1> [0] 1= Use clovr defaults for file paths 
		--diag=				<0|1> [0] 1= Use diag node defaults for file paths 
		--fs=				<0|1> [0] 1= Use filesystem defaults for file paths 
		--help				Help Basic Info
		--help_full 			Help Full Info\n";
}


if(!$options{input}){die "Error: Please give an input.bam with --input=<FILE>. Try again or use --help.\n";}
if(!$options{output_dir}){print "It is HIGHLY recommended you STOP, restart, and use a --output_dir=<some/where/>.\n";sleep 60;}

## Setup Default paths for references and bins:
my $lgtseek = LGTSeek->new2({
	options => \%options,
	});

my $inputs = setup_input(\%options);
foreach my $input (@$inputs){
	if($lgtseek->{decrypt}==1 && !$lgtseek->{url}){die "Error: Must give a --url to use --decrypt.\n";}
	my ($name,$path,$suf)=fileparse($lgtseek->{input},('.gpg.bam','_prelim.bam','.bam'));
	chomp $name;
	if(!$lgtseek->{output_dir}){$lgtseek->{output_dir}="$path/lgtseq/";}
	if($lgtseek->{subdirs}==1){$lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}"); $lgtseek->{output_dir} = "$options{output_dir}/"."$name/";}

	## Qsub the job instead of running it
	if($lgtseek->{Qsub}==1){
		my $cmd = "/home/ksieber/scripts/lgt_seq.pl";
		foreach my $key (keys %options){
			next if ($key=~/Qsub/);
			if($options{$key}){$cmd = $cmd." --$key=$options{$key}"};
		}
		Qsub2({
			cmd => $cmd,
			threads => "$lgtseek->{threads}",
			wd => "$lgtseek->{output_dir}",
			project => "$lgtseek->{project}",
			});
		last;
	}


	my $print = "lgt_seq.pl";
	foreach my $key (keys %options){if($options{$key}){$print = "$print"." \-\-$key=$options{$key}";}}
	print STDERR "\n$print\n";

	print STDERR "\n++++++++++++++++++++++++++++\n";
	print STDERR "++++++++  LGT-SEQ  +++++++++\n";
	print STDERR "++++++++++++++++++++++++++++\n\n";

	# Get input ready for processing with decryption and/or prelim filtering
	print STDERR "=========PREP-INPUT========\n";
	my $input_bam;
	if($lgtseek->{decrypt}==1){
		$input_bam = $lgtseek->decrypt({
			input => $lgtseek->{input},
			url => $lgtseek->{url},
			output_dir => $lgtseek->{output_dir}
			});
	} elsif ($lgtseek->{prelim_filter}==1) {
		my $unfiltered_bam;
		if ($lgtseek->{decrypt}==1){ $unfiltered_bam = $input_bam; }
		else { $unfiltered_bam = $lgtseek->{input}; }
		my $bams_array_ref = $lgtseek->prelim_filter({
			input_bam => $unfiltered_bam,
			output_dir => "$lgtseek->{output_dir}/prelim_filter/",
			overwrite => 0,
			});
		$input_bam = ${$bams_array_ref}[0];
	} else { $input_bam = $lgtseek->{input}; }
	print STDERR "=========INPUT-READY=======\n";

	# Align to the donors.
	print STDERR "========RUNBWA-DONOR========\n";
	my $donor_bams = $lgtseek->runBWA({
		input_bam => $input_bam,
		output_bam => 1,
		threads => $lgtseek->{threads},
		output_dir => "$lgtseek->{output_dir}/donor_alignments/",
		reference_list => $lgtseek->{split_bac_list},
		overwrite => 0,   
		cleanup_sai => 1,
		});

	time_check();

	# Align to the hosts.
	print STDERR "========RUNBWA-HOST========\n";
	my $host_bams = $lgtseek->runBWA({
		input_bam => $input_bam,
		output_bam => 1,
		threads => $lgtseek->{threads},
		output_dir => "$lgtseek->{output_dir}/host_alignments/",
		reference => $lgtseek->{hg19_ref},
		overwrite => 0, 
		cleanup_sai => 1,   
		});
	time_check();

	# Postprocess the results
	print STDERR "=======POSTPROCESS=======\n";
	my $pp_data = $lgtseek->bwaPostProcess({
		donor_bams => $donor_bams,
		host_bams => $host_bams,
		output_prefix => $name,
		overwrite => 0,   
		});
	time_check();
	# Clean up output we don't need anymore; This is IMPORTANT on nodes.
	print STDERR "RM: $lgtseek->{output_dir}/host_alignments/\n";
	print STDERR $lgtseek->_run_cmd("rm -rf $lgtseek->{output_dir}/host_alignments/");
	print STDERR "RM: $lgtseek->{output_dir}/donor_alignments/\n";
	print STDERR $lgtseek->_run_cmd("rm -rf $lgtseek->{output_dir}/donor_alignments/");

	# Create file with number of counts
	my @header = ('run_id');
	my @vals = ($name);
	open OUT, ">$lgtseek->{output_dir}/$name\_post_processing.tab" or die;
	map {
		push(@header,$_);
		my $foo = $pp_data->{counts}->{$_} ? $pp_data->{counts}->{$_} : 0;
		push(@vals,$foo);
	} ('total','host','no_map','all_map','single_map','integration_site_host','integration_site_donor','microbiome','lgt');
	&print_tab("$lgtseek->{output_dir}/$name\_post_processing.tab",\@header,\@vals);


	## Check to make sure we found LGT.
	print STDERR "========LGT========\n";
	if($lgtseek->empty_chk({input => $pp_data->{files}->{lgt_donor}})==1){
		print STDERR "No LGT in: $pp_data->{files}->{lgt_donor}\. Skipping LGT LCA calculation and blast validation.\n";
	} else {
		# Prinseq filter the putative lgts
		print STDERR "=======LGT-PRINSEQ=======\n";
		my $filtered_bam = $lgtseek->prinseqFilterBam({
			output_dir => "$lgtseek->{output_dir}/lgt_prinseq_filtering",
			input_bam => $pp_data->{files}->{lgt_host}
			});

		# Add filtered count to counts.
		push(@header,'lgt_pass_prinseq');
		push(@vals,$filtered_bam->{count});
		&print_tab("$lgtseek->{output_dir}/$name\_post_processing.tab",\@header,\@vals);

		# Calculate BWA LCA's for LGTs
		print STDERR "======LGT-BWA-LCA======\n";
		print STDERR `mkdir -p $lgtseek->{output_dir}/lgt_lca-bwa`;
		$lgtseek->runBWA({
			input_bam => "$lgtseek->{output_dir}\/$name\_lgt_host_filtered.bam",
			output_dir => "$lgtseek->{output_dir}\/lgt_lca-bwa/",
			out_file => "$lgtseek->{output_dir}\/lgt_lca-bwa\/$name\_lgt_lca-bwa.txt",
			reference_list => "$lgtseek->{refseq_list}",
			threads => $lgtseek->{threads},
			run_lca => 1,
			overwrite => 0,
			cleanup_sai => 1,
			});
		# Make LGT Fasta for Blast validation
		my $lgt_fasta = $lgtseek->sam2Fasta({
			input => "$lgtseek->{output_dir}\/$name\_lgt_host_filtered.bam",
			output_dir => "$lgtseek->{output_dir}/blast_validation/"
			});
		time_check();

		# Blast & get best hits
		print STDERR "=======BESTBLAST2=======\n";
		my $best_blasts = $lgtseek->bestBlast2({
			db => $lgtseek->{path_to_blastdb},
			lineage1 => $lgtseek->{donor_lineage},
			lineage2 => $lgtseek->{host_lineage},
			fasta => $lgt_fasta,
			output_dir => "$lgtseek->{output_dir}/blast_validation/"
			});
		time_check();

		# Now run lgtfinder
		print STDERR "========LGTFINDER=========\n";
		my $valid_lgts = $lgtseek->runLgtFinder({
			lineage1 => $lgtseek->{donor_lineage},
			lineage2 => $lgtseek->{host_lineage},
			input_file_list => $best_blasts->{list_file},
			output_prefix => "$name",
			output_dir => "$lgtseek->{output_dir}/lgt_finder/",
			});
		time_check();

		# Run blast and keep raw output ?
		my $blast_ret = $lgtseek->_run_cmd("blastall -p blastn -e 10e-5 -T F -d $lgtseek->{path_to_blastdb} -i $lgt_fasta > $lgtseek->{output_dir}/blast_validation/$name\_blast.raw");
		
		## Cleanup Intermediate Files
		print STDERR "RM: $lgtseek->{output_dir}/lgt_prinseq_filtering/\n";
		$lgtseek->_run_cmd("rm -rf $lgtseek->{output_dir}/lgt_prinseq_filtering/");
		
		time_check();
	}

	# Calculate BWA LCA's for Microbiome Reads
	print STDERR "======Microbiome======\n";
	## Check to make sure we found Microbiome Reads. If no microbiome reads skip this step. 
	if($lgtseek->empty_chk({input => "$lgtseek->{output_dir}\/$name\_microbiome.bam"})==1){
		print STDERR "No Microbiome reads in: $lgtseek->{output_dir}\/$name\_microbiome.bam. Skipping microbiome LCA calculation.\n";
	} else {
		# Prinseq filter the putative lgts
		print STDERR "=====Microbiome-PRINSEQ=====\n";
		my $filtered_bam = $lgtseek->prinseqFilterBam({
			output_dir => "$lgtseek->{output_dir}/microbiome_prinseq_filtering",
			input_bam => "$lgtseek->{output_dir}\/$name\_microbiome.bam",
			});
		time_check();

		# Add filtered count to counts.
		push(@header,'microbiome_pass_prinseq');
		push(@vals,$filtered_bam->{count});
		&print_tab("$lgtseek->{output_dir}/$name\_post_processing.tab",\@header,\@vals);

		print STDERR "====Microbiome-BWA-LCA====\n";
		$lgtseek->runBWA({
			input_bam => "$lgtseek->{output_dir}\/$name\_microbiome_filtered.bam",
			output_dir => "$lgtseek->{output_dir}\/microbiome_lca-bwa/",
			out_file => "$lgtseek->{output_dir}\/microbiome_lca-bwa\/$name\_microbiome_lca-bwa.txt",
			reference_list => "$lgtseek->{refseq_list}",
			threads => $lgtseek->{threads},
			run_lca => 1,
			overwrite => 0,
			cleanup_sai => 1,
			});

		## Cleanup Intermediate Files
		print STDERR "RM: $lgtseek->{output_dir}/microbiome_prinseq_filtering/\n";
		$lgtseek->_run_cmd("rm -rf $lgtseek->{output_dir}/microbiome_prinseq_filtering/");

		time_check(); 
	}

	# Calculate coverage of LGT on human side.
	if($lgtseek->{lgt_coverage}==1){
		print STDERR "======Calculating Coverage of Hg19 LGT======\n";
		$lgtseek->mpileup({
			input => "$lgtseek->{output_dir}\/$name\_lgt_host_filtered.bam",
			output_dir => $lgtseek->{output_dir},
			ref => $lgtseek->{hg19_ref},
			cleanup => 1,
			overwrite => 0,
			});
		time_check();
	}

	time_check();
	print STDERR "======Completed lgt_seq.pl on $lgtseek->{input}======\n";
}

sub print_tab {
	my ($file,$header,$vals) = @_;
	open OUT, ">$file" or die "Couldn't open $file\n";
	print OUT join("\t",@$header);
	print OUT "\n";
	print OUT join("\t",@$vals);
	print OUT "\n";
}

sub time_check {
	my $elapsed = runtime(); 
	print STDERR "Time Since Start: $elapsed\n";
}



__END__
