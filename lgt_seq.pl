#!/usr/bin/perl

=head1 NAME

lgt_seq.pl

=head1 SYNOPSIS

Search an bam for bacterial-human LGT.

=head1 DESCRIPTION


=head1 AUTHORS - Karsten Sieber

e-mail: Karsten.sieber@gmail.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with "_"

=cut
use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
		'input|i=s', # Comma separated list of files
		'input_list|I=s',
		'Qsub|Q=i',
		'sub_mem|mf=s',
		'sub_name=s',
		'excl=i',
		'no_gal=i',
		'hostname=s',
		'decrypt=i',
		'url=s',
		'prelim_filter=i',
		'name_sort_input=i',
		'keep_softclip=i',
		'split_bam=i',
		'seqs_per_file=i',
		'aln1_human=i',
		'aln2_human=i',
		'split_bac_list=s',
		'hg19_ref=s',
		'refseq_list=s',
		'output_dir|o=s',
		'subdirs=i',
		'lgt_coverage=i',
		'max_overlap=i',
		'min_length=i',
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
		'clovr=i',
		'diag=i',
		'fs=i',
		'verbose|V=i',
		'print_hostname|ph=i',
		'config_file=s',
		'help|h',
		'help_full|H',
		'workflow_help'
		);
use print_call;
print_hostname(\%options);									## This is useful for trouble shooting grid nodes that might be missing modules for LGTSeek etc. 
use lib '/local/projects-t3/HLGT/scripts/lgtseek/lib/';      ### May need to change this depending on where the script is being run
use LGTSeek;
use Time::SoFar;
use run_cmd;
use setup_input;
use Carp;
$Carp::MaxArgLen = 0;
use File::Basename;

## Check if the user needs help information
if($options{help}){&help;}								## @ end of script
if($options{help_full}){&help_full;}					## @ end of script
if($options{workflow_help}){&workflow_help;}			## @ end of script

## Check we have the necessary inputs
if(!$options{input} && !$options{input_list}){confess "Error: Please give an input.bam with --input=<FILE> or --input_list=<LIST>. Try again or use --help_full.\n";}
if(!$options{output_dir}){print "It is HIGHLY recommended you STOP, restart, and use a --output_dir=</some/path/>.\n";sleep 60;}
if(!$options{sub_name}){$options{sub_name}="lgtseq";}

## Qsub the job instead of running it
if($options{Qsub}){Qsub3(\%options)};

## Print the script call
print_call(\%options);
## Initialize LGTSeek.pm
my $lgtseek = LGTSeek->new2(\%options);
## Setup array ref of inputs
my $inputs = setup_input(\%options);

foreach my $input (@$inputs){
	my ($name,$path,$suf) = fileparse($input,$lgtseek->{suffix_regex});
	chomp $name;
	## Setup output directory
	if(!$lgtseek->{output_dir}){$lgtseek->{output_dir}="$path/lgtseq/";}
	if($lgtseek->{subdirs}){
		$lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}"); 
		$lgtseek->{output_dir} = "$options{output_dir}/"."$name/"; 
		$lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}");
	}

	# Primary aln to Human.
	if($lgtseek->{aln1_human}){ print STDERR "======== RUNBWA-Human1 ========\n";
		my $human_bam1;
		## Map input bam @ hg19
		if($input =~/\.bam$/){
			$human_bam1 = $lgtseek->runBWA({								## &runBWA returns an array
				input_bam => $input,
				output_bam => 1,
				threads => $lgtseek->{threads},
				output_dir => "$lgtseek->{output_dir}/human_aln1/",
				reference => $lgtseek->{hg19_ref},
				overwrite => $lgtseek->{overwrite},   
				cleanup_sai => 1,
				});
		}
		## Map fastqs @ hg19
		elsif($input=~/.fq$/ || $input=~/.fastq$/ || $input=~/.fastq.gz$/){
			my ($in1,$in2) = split(/,/,$input);								## split input if needed
			my ($input_base,$input_dir,$input_suffix)=fileparse($in1,@{$lgtseek->{fastq_suffix_list}});
			$human_bam1 = $lgtseek->runBWA({								## &runBWA returns an array
				input_dir => $input_dir,
				input_base => $input_base,
				output_bam => 1,
				threads => $lgtseek->{threads},
				output_dir => "$lgtseek->{output_dir}/human_aln1/",
				reference => $lgtseek->{hg19_ref},
				overwrite => $lgtseek->{overwrite},   
				cleanup_sai => 1,
				});
		}
		$input = @$human_bam1[0];
		time_check();
	} 

	# Prelim_filter
	if($lgtseek->{prelim_filter} || $lgtseek->{name_sort_input}){ print STDERR "======== Prelim_filter ========\n";
		my $prelim_filtered_bam = $lgtseek->prelim_filter({				## &prelim_filter returns an array
			input_bam => $input,
			output_dir => "$lgtseek->{output_dir}/prelim_filter/",		## vv lgtseek &prelim_filter defaults vv
			name_sort_input => $lgtseek->{name_sort_input},				## Default = 1
			sort_mem => $lgtseek->{sort_mem},							## Default = 10G lgtseek default. lgt_prep overides to 40G. 
			threads => $lgtseek->{threads},								## Default = 1
			split_bam => "0",											## Default = 1
			keep_softclip => $lgtseek->{keep_softclip},					## Default = 1, better to split with lgt_seq_prelim
			overwrite => $lgtseek->{overwrite},							## Default = 0
			});															## ^^								^^
		$input = @$prelim_filtered_bam[0];								## This works because we are not splitting bams
	} 

	# Secondary aln human
	if($lgtseek->{aln2_human}){ print STDERR "======== RUNBWA-Human2 ========\n";
		my $human_bam2 = $lgtseek->runBWA({								## &runBWA returns an array
			input_bam => $input,
			output_bam => 1,
			threads => $lgtseek->{threads},
			output_dir => "$lgtseek->{output_dir}/human_aln2/",
			reference => $lgtseek->{hg19_ref},
			overwrite => $lgtseek->{overwrite},   
			cleanup_sai => 1,
			});
		$input = @$human_bam2[0];
		time_check();
	}

	# Align to the Bacteria.
	print STDERR "======== RUNBWA-Bacteria ========\n";
	my $bacterial_bams = $lgtseek->runBWA({								## &runBWA returns an array
		input_bam => $input,
		output_bam => 1,
		threads => $lgtseek->{threads},
		output_dir => "$lgtseek->{output_dir}/bacterial_alignments/",
		reference_list =>  $lgtseek->{split_bac_list},
		overwrite => $lgtseek->{overwrite}, 
		cleanup_sai => 1,   
		});
	time_check();

	# Postprocess the BWA mappings to human and bacteria
	print STDERR "======= POSTPROCESS =======\n";
	my $host_bams = ["$input"]; 										## Put human input bam into array ref for post-proc.
	my $pp_data = $lgtseek->bwaPostProcess({
		donor_bams => $bacterial_bams,
		host_bams => $host_bams,
		output_prefix => $name,
		overwrite => $lgtseek->{overwrite},   
		});
	time_check();

	# Create file with number of counts
	my @header = ('run_id');
	my @vals = ($name);
	open OUT, ">$lgtseek->{output_dir}/$name\_post_processing.tab" or die;
	map {
		push(@header,$_);
		my $foo = $pp_data->{counts}->{$_} ? $pp_data->{counts}->{$_} : 0;
		push(@vals,$foo);
	} ('total','host','no_map','all_map','single_map','integration_site_human','integration_site_bac','microbiome','lgt');
	&print_tab("$lgtseek->{output_dir}/$name\_post_processing.tab",\@header,\@vals);

	# Clean up output we don't need anymore; This is IMPORTANT on nodes. Not so much on filesystem.
	$lgtseek->_run_cmd("rm -rf $lgtseek->{output_dir}/human_aln1/");
	$lgtseek->_run_cmd("rm -rf $lgtseek->{output_dir}/human_aln2/");
	$lgtseek->_run_cmd("rm -rf $lgtseek->{output_dir}/prelim_filter/");
	$lgtseek->_run_cmd("rm -rf $lgtseek->{output_dir}/bacterial_alignments/");
	
	##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	## LGT Analysis.
	##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	## Check to make sure we found LGT.
	print STDERR "======== LGT ========\n";
	if($lgtseek->empty_chk({input => $pp_data->{files}->{lgt_donor}})){
		print STDERR "No LGT in: $pp_data->{files}->{lgt_donor}\. Skipping LGT LCA calculation and blast validation.\n";
	} else {
		# Prinseq filter the putative lgts
		# This is for manual curation later. Use all potential LGT for further analysis.
		print STDERR "======= LGT-PRINSEQ =======\n";
		my $filtered_bam = $lgtseek->prinseqFilterBam({
			output_dir => $lgtseek->{output_dir},
			input_bam => $pp_data->{files}->{lgt_host}
			});

		# Add filtered count to counts.
		push(@header,'lgt_pass_prinseq');
		push(@vals,$filtered_bam->{count});
		&print_tab("$lgtseek->{output_dir}/$name\_post_processing.tab",\@header,\@vals);

		# Calculate BWA LCA's for LGTs
		print STDERR "====== LGT-BWA-LCA ======\n";
		$lgtseek->runBWA({
			input_bam => $filtered_bam->{bam},
			output_dir => "$lgtseek->{output_dir}\/lgt_bwa-lca/",
			out_file => "$lgtseek->{output_dir}\/lgt_bwa-lca\/$name\_lgt_lca-bwa.txt",
			reference_list => $lgtseek->{refseq_list},
			threads => $lgtseek->{threads},
			run_lca => 1,
			overwrite => $lgtseek->{overwrite},
			cleanup_sai => 1,
			});

		# Blast & get best hits
		print STDERR "========== LGT-BESTBLAST2 ==========\n";
		my $best_blasts = $lgtseek->bestBlast2({
			db => $lgtseek->{path_to_blastdb},
			lineage1 => $lgtseek->{donor_lineage},
			lineage2 => $lgtseek->{host_lineage},
			bam => $filtered_bam->{bam},
			output_dir => "$lgtseek->{output_dir}/blast_validation/"
			});
		time_check();

		# LGTFinder
		print STDERR "======== LGTFINDER =========\n";
		my $valid_lgts = $lgtseek->runLgtFinder({
			lineage1 => $lgtseek->{donor_lineage},
			lineage2 => $lgtseek->{host_lineage},
			input_file_list => $best_blasts->{list_file},
			output_prefix => $name,
			max_overlap => $lgtseek->{max_overlap},
			min_length => $lgtseek->{min_length},
			output_dir => "$lgtseek->{output_dir}/lgt_finder/",
			});
	
		## Create a new bam from blast validation & lgtfinder results && add #'s to post_processing.tab		
		my $blast_validated_lgt_bam = $lgtseek->validated_bam({ 
			input => $filtered_bam->{bam}, 
			by_clone => $valid_lgts->{by_clone}, 
			output => "$lgtseek->{output_dir}/$name\_lgt_host_filtered_validated.bam" 
			});

		## Add # validated to post_processing.tab
		push(@header,'lgt_valid_blast');
		push(@vals,"$blast_validated_lgt_bam->{count}");
		&print_tab("$lgtseek->{output_dir}/$name\_post_processing.tab",\@header,\@vals);
		time_check();

		# Run blast and keep raw output
		#my $blast_ret = $lgtseek->_run_cmd("blastall -p blastn -e 10e-5 -T F -d $lgtseek->{path_to_blastdb} -i $lgt_fasta > $lgtseek->{output_dir}/blast_validation/$name\_blast.raw");	## FIX KBS 01.06.14
		#time_check();
	}

	##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	## Microbiome Analysis
	##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	## Calculate BWA LCA's for Microbiome Reads
	print STDERR "====== Microbiome ======\n";
	## Check to make sure we found Microbiome Reads. If no microbiome reads skip this step. 
	if($lgtseek->empty_chk({input => "$lgtseek->{output_dir}\/$name\_microbiome.bam"})){
		print STDERR "No Microbiome reads in: $lgtseek->{output_dir}\/$name\_microbiome.bam. Skipping microbiome LCA calculation.\n";
	} else {
		# Prinseq filter the putative microbiome reads
		# Use only prinseq quality reads for microbiome analysis
		print STDERR "===== Microbiome-PRINSEQ =====\n";
		my $filtered_bam = $lgtseek->prinseqFilterBam({
			input_bam => $pp_data->{files}->{microbiome_donor},															## KBS 01.08.14 "$lgtseek->{output_dir}\/$name\_microbiome.bam",
			output_dir => $lgtseek->{output_dir},
		});
		time_check();

		# Add filtered count to counts.
		push(@header,'microbiome_pass_prinseq');
		push(@vals,$filtered_bam->{count});
		&print_tab("$lgtseek->{output_dir}/$name\_post_processing.tab",\@header,\@vals);

		print STDERR "==== Microbiome-BWA-LCA ====\n";
		$lgtseek->runBWA({
			input_bam => $filtered_bam->{bam},
			output_dir => "$lgtseek->{output_dir}\/microbiome_bwa-lca/",
			out_file => "$lgtseek->{output_dir}\/microbiome_bwa-lca/$name\_microbiome_lca-bwa.txt",					 ## KBS 01.07.14
			reference_list => $lgtseek->{refseq_list},
			threads => $lgtseek->{threads},
			run_lca => 1,
			overwrite => $lgtseek->{overwrite},
			cleanup_sai => 1,
		});

		## Cleanup Intermediate Files
		print STDERR "RM: $lgtseek->{output_dir}/microbiome_prinseq_filtering/\n";
		$lgtseek->_run_cmd("rm -rf $lgtseek->{output_dir}/microbiome_prinseq_filtering/");
		time_check(); 
	}

	##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	## Calculate coverage of LGT on human side. (Not apropriate most the time)
	##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if($lgtseek->{lgt_coverage}==1){
		print STDERR "====== Calculating Coverage of Hg19 LGT ======\n";
		$lgtseek->mpileup({
			input => "$lgtseek->{output_dir}\/$name\_lgt_host_filtered.bam",
			output_dir => $lgtseek->{output_dir},
			ref => $lgtseek->{hg19_ref},
			cleanup => 1,
			overwrite => $lgtseek->{overwrite},
		});
		time_check();
	}

	time_check();
	print_complete(\%options);
}




## Subroutines
sub print_tab {
	my ($file,$header,$vals) = @_;
	open OUT, ">$file" or $lgtseek->fail("Couldn't open $file\n");
	print OUT join("\t",@$header);
	print OUT "\n";
	print OUT join("\t",@$vals);
	print OUT "\n";
}

sub time_check {
	my $elapsed = runtime(); 
	print STDERR "Time Since Start: $elapsed\n";
}

sub help {
	die "Help: This script will identify bacterial human LGT.
		--input=			<Input> Accepts bams, fastq, and fastq.gz. With fatsq's only use 1 of the pair for input. (ie: foo_1.fastq.gz)
		--output_dir=			Directory for all output. Will be created if it doesn't exist. 
		--help_full 			Full help info on options.
		--workflow_help			Help setting up an efficient lgtseq workflow with the optional steps.\n";
}

sub help_full {
	die "Help: This script takes a bam and identifies bacterial human LGT.
	     _____
	____/Input\\_________________________________________________________________________________
		--input|i=			<Input BAM or fastq> Needs to be name sorted. Also use --prelim_filter & --name_sort_input for position sorted inputs.
		--input_list|I=			<List of inputs> 1 per line.
	     ______
	____/Output\\________________________________________________________________________________
		--output_dir|o=			Directory for all output. Will be created if it doesn't exist. 
	  	  --subdirs=			<0|1> [0] 1= Make a sub-directory in output_dir based on input name
	     _______________________
	____/Primary Human Alignment\\_______________________________________________________________
		--aln1_human=			<0|1> [0] 1= Primary aln to hg19.
	     ___________________
	____/Prelimary Filtering\___________________________________________________________________
		--prelim_filter=		<0|1> [0] 1= Filter a human mapped bam, keeping potential LGT & Microbiome reads.
	  	  --name_sort_input= 		<0|1> [0] 1= Resort the input bam by read names.
		  --keep_softclip=		<0|1> [1] 1= Keep soft clipped reads >=24 bp (Reads potentially on LGT)
	  	  --sort_mem=			[5G] Mem per thread to sort with. Careful this corresponds with --threads. 
	     _________________________
	____/Secondary Human Alignment\\______________________________________________________________
		--aln2_human=			<0|1> [0] 1= Secondary aln to hg19. Useful for mapping to standarized hg19 after prelim_filter.
	     __________________
	____/Submit to SGE-grid\\_____________________________________________________________________
		--Qsub|Q=			<0|1> [0] 1= qsub the job to the grid.
		  --threads=|t			[1] # of CPU's to use for multithreading BWA sampe
		  --sub_mem|mf=			[6G] Min mem to qsub for on grid
		  --sub_name=			[lgtseq] 
		  --project=			Grid project to use. 
		  --no_gal=			<0|1> [1] Avoid galactus node. 
		  --excl=			<0|1> [1] Run exclusively. (Not suggested)
		  --print_hostname|ph= 		<0|1> [1] Print hostname. Defaults to 1 if verbose=1.
	     __________
	____/References\\_____________________________________________________________________________
		--hg19_ref=			Path to hg19 reference
		--split_bac_list=		Path to the list of split bacterial references (A2D, E2P, R2Z)
		--refseq_list=			Path to all bacterial references in refseq.  
	     _______________
	____/Bin Directories\\________________________________________________________________________
		--bin_dir=			[/local/projects-t3/HLGT/scripts/lgtseek/bin/]
		--taxon_host=			[mongotest1-lx.igs.umaryland.edu:10001]
		--taxon_dir=			[/local/db/repository/ncbi/blast/20120414_001321/taxonomy/]
		--taxon_idx_dir=		[/local/projects-t3/HLGT/idx_dir/20120414]
		--path_to_blastdb=		[/local/db/repository/ncbi/blast/20120414_001321/nt/nt]
	     ________________
	____/Help Information\\_______________________________________________________________________
		--verbose|V			<0|1> [1] 1= Verbose reporting of progress. 0 =Turns off reports. 
		--help|h			Help Basic Info
		--help_full|H			Help Full Info
		--workflow_help			Examples how to use optional portions of lgt_seq to increase efficiency.
	_____________________________________________________________________________________________\n";
}

sub workflow_help {
	die 
"----------------------------------------------------------------------------------------
Workflow of modules:
	input -> Primary aln hg19 -> Prelim_filter -> Secondary aln hg19 -> map_bacteria -> LGTFINDER
	Options: --aln1_human=0|1 --prelim_filter=0|1  --aln2_human=0|1	Always 1     Always 1
	Ex1 - Input bacteria mapped:	 		--aln1_human=1 --prelim_filter=1 --aln2_human=0
	Ex2 - Input human mapped 2 non-hg19: 		--aln1_human=0 --prelim_filter=1 --aln2_human=1
	Ex3 - Input position sorted	         	--aln1_human=0 --prelim_filter=1 --aln2_human=0
----------------------------------------------------------------------------------------\n";
}

