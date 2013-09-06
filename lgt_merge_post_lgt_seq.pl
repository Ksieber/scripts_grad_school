#!/usr/bin/perl
use warnings;
use strict;
use lib qw(/local/projects-t3/HLGT/scripts/lgtseek/lib/ /opt/lgtseek/lib/);      ### May need to change this depending on where the script is being run
use LGTSeek;
use run_cmd;
use setup_input;
use Time::SoFar;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
	'input=s',
	'input_list=s',
	'output_dir=s',
	'subdirs=s',
	'help',
	);

if($options{help}){
	die "Help:
	This script will merge the output files of lgt_seq.pl.
	----------------------------------------------------------------------------------------
	--input=		Directory to merge files from. Subdirectories in this directory should contain lgt_seq.pl outputs.
	--input_list=	List of input directories.
	----------------------------------------------------------------------------------------
	--output_dir=	Directory for output.
	  --subdirs=		Make a new directory within output_dir for each input directory.
	----------------------------------------------------------------------------------------
	--help
	----------------------------------------------------------------------------------------\n";
}

if(!$options{input} && !$options{input_list}){die "Must give an input. Use --input or --input_list\n";}
if(!$options{output_dir}){die "Must use --output_dir=\n";}
my $input = setup_input(\%options);
run_cmd("mkdir -p $options{output_dir}");

my @bam_suffix_to_merge = (
		'microbiome.bam',
		'lgt_host.bam',
		'integration_site_donor_host.bam',
		'lgt_host_filtered.bam',
		'integration_site_donor_donor.bam',
		'microbiome_filtered.bam'
		);

my @txt_suffix_to_merge = (
		'by_clone.txt',
		'by_trace.txt',
		'post_processing.tab',
		'lgt_host_lineage1.out',
		'lgt_host_lineage2.out',
		'microbiome_lca-bwa_independent_lca.txt',
		'microbiome_lca-bwa.txt',
		'lgt_lca-bwa_independent_lca.txt',
		'lgt_lca-bwa.txt',
		);

foreach my $dir (@$input){
	$dir =~ /\/(\w+)\/$/;
	my $name = $1;
	my $output_dir;
	if($options{subdirs}){
		run_cmd("mkdir -p $options{output_dir}/$name");
		$output_dir="$options{output_dir}/$name";
	} else {
		$output_dir="$options{output_dir}";
	}

	foreach my $txt_suffix (@txt_suffix_to_merge){
		chomp(my @list_to_merge = `find $dir -name '*$txt_suffix'`);
		my $output = "$options{output_dir}/$name\_$txt_suffix";
		foreach my $file (@list_to_merge){
			run_cmd("cat $file >> $output");
		}
	}

	foreach my $bam_suffix (@bam_suffix_to_merge){
		chomp(my @list_to_merge = `find $dir -name '*$bam_suffix'`);
		my $output = "$options{output_dir}\/$name\_$bam_suffix";
		my $header = run_cmd("samtools view -H $list_to_merge[0]");
		open(my $out, "| samtools view -S - -bo $output") or die "Can not open output: $output\n";		
		print $out "$header\n";
		foreach my $bam (@list_to_merge){
			open(my $in, "-|","samtools view $bam") or die "Can not open input: $bam\n";
			while(<$in>){
				print $out "$_";
			}
			close $in or die "Can not close input: $input because: $!\n";
		}
		close $out or die "Can't close output: $output because: $!\n";
	}
}
