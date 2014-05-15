#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
## Implement gff -> draw better "reference contig"
## Implement bio::Fetch to bio:seq bam files?
## Implement bio::seq obj -> draw multiple reads?
use warnings;
use strict;
use lib "/local/projects-t3/HLGT/scripts/lgtseek/lib/";      ### May need to change this depending on where the script is being run
use LGTSeek;
use Bio::Graphics;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::DB::Sam;
use File::Basename;
use run_cmd;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
		'input|i=s',
		'region=s',
		'sort_bam=s',
		'scale=s',
		'region=s',
		'image_length=s',
		'image_width=s',
		'pad_scale=s',
		'output_prefix=s',
		'output_dir|o=s',
		'help',
		'help_full'
		) or die "Error: Unrecognized command line option. Please try again.\n";
if($options{help}){&help;}
if($options{help_full}){&help_full;}

my $lgtseek = LGTSeek->new2(\%options);
#############################################################
## Setup & check input
if(!$options{input} && !$ARGV[0]){die "Must give an input: --input=<BAM> or ARGV[0]\n";}
my $input = defined $options{input} ? $options{input} : $ARGV[0];
my ($fn,$path,$suff)=fileparse($input,('.srt.bam','.bam'));
if($lgtseek->empty_chk({input => $input})==1){die "Error: The input: $options{input} is empty.\n";}
#############################################################
## Setup Output
my $output_dir = defined $options{output_dir} ? "$options{output_dir}" : "$path";
my $output_prefix =  defined $options{output_prefix} ? "$options{output_prefix}" : "$fn";
my $out = "$output_dir/$output_prefix\.png";
#############################################################
## Setup Defaults
my $sort_bam = defined $options{sort_bam} ? $options{sort_bam} : "0";
my $region = defined $options{region} ? $options{region} : undef;							## PROBLEM
my $image_width = defined $options{image_width} ? $options{image_width} : "1000";
my $image_length = defined $options{image_length} ? $options{image_length} : "500";
my $pad_scale = defined $options{pad_scale} ? $options{pad_scale} : "0";
#############################################################
## Sort bam
if($sort_bam==1){
	run_cmd("samtools sort $input $output_dir/$fn\_bam2png_tmp_srt");
	run_cmd("samtools index $fn\_bam2png_tmp_srt.bam $output_dir/$fn\_bam2png_tmp_srt.bai");
	$input = "$output_dir/$fn\_bam2png_tmp_srt.bam";	
}
#############################################################
## Calculate scale size 
my $first_position = defined $options{region} ? run_cmd("samtools view -hu $input $options{region} | samtools mpileup -A - | head -n 1 | cut -f2") : run_cmd("samtools mpileup -A $input | head -n 1 | cut -f2");
my $last_position  = defined $options{region} ? run_cmd("samtools view -hu $input $options{region} | samtools mpileup -A - | tail -n 1 | cut -f2") : run_cmd("samtools mpileup -A $input | tail -n 1 | cut -f2");
my $scale_size = defined $options{scale} ? $options{scale} : (($first_position - $pad_scale )+($last_position + $pad_scale));		## = Total + padding on each side
#############################################################
## Setup the graphics panels and tracks.
## Create the Panel
my $panel = Bio::Graphics::Panel->new(
					-length    => $image_length,
					-width     => $image_width,
					-pad_left  => 10,
					-pad_right => 10,
					);
## Scale at the top of the img
my $scale = Bio::SeqFeature::Generic->new(
	-start => 1,								## Need to make this the start position
	-end   => $scale_size,						## Need to make this the end position
	);
$panel->add_track($scale,
	-glyph   => 'arrow',
	-tick    => 2,
	-fgcolor => 'black',
	-double  => 1,
	);
## Track for bam data 
my $track = $panel->add_track(
					-glyph     => 'transcript2',
					-label     => 1,
					-bgcolor   => 'green',
					);

#############################################################
my $bam = Bio::DB::Sam->new(
	-bam => $input,
	-fasta => $options{ref},
	-expand_flags => 1,
	-autoindex => 1);

## Add feature by read_list and region later
my @pairs = $bam->features(-type => 'read_pair');
for my $pair (@pairs){
	my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
	$track->add_feature([$first_mate,$second_mate]);
}

# $track->add_feature(arrow => \@pairs);

#############################################################
my $OFH;
if($options{output_dir} || $options{output_prefix}){
	open($OFH,">","$out") or die "Can't output the output: $out because: $!\n";
} else {
	open ($OFH, " | display - ") or die "Can't output to display.\n";
}
print $OFH $panel->png;
close $OFH;



sub help {
	die 
"	----------------------------------------------------------------------------------------
	Help: This script will take a bam and use Bio::Graphics to create a .png of the reads.
		----------------------------------------------------------------------------------------
		--input=			Position Sorted bam. Should be very small b/c memory intense. 
		--help_full			Complete help information.
		----------------------------------------------------------------------------------------
";
}

sub help_full {
	die	
"	----------------------------------------------------------------------------------------
	Help: This script will take a bam and use Bio::Graphics to create a .png of the reads.
		The bam either needs to cover a small region or --region should be used. 
		ex: bam2png.pl foo.bam | display - 
		----------------------------------------------------------------------------------------
		## ARGV[0]				Position sorted bam. Should be very small b/c memory intense. 
		--input=			Position Sorted bam. Should be very small b/c memory intense. 
		--sort_bam=			<0|1> [0] 1= Sort the input bam on position.
		----------------------------------------------------------------------------------------
		--region=			Region of the sorted bam to draw. ex: \"chr19:100-200\" (Must follow this format)
		--scale=			Size of the scale to draw on the top of the image. [region size]
		--image_length=			Length of the image. [1000]
		--image_width= 			Width of the image. [500]
		--pad_scale=			Number of bp to pad eachside with. [0] ## SHOULD MAKE IT GO NEGATIVE FROM START 
		----------------------------------------------------------------------------------------
		Default output =>	 	STDOUT. Display - || > foo.png
		--output_prefix=		[filename] Prefix for the output file. Optional.
		--output_dir=			[input_dir] Directory to put the output. Optional.
		----------------------------------------------------------------------------------------
";
}

__END__
