#!/usr/bin/perl
use warnings;
use strict;
use lib "/local/projects-t3/HLGT/scripts/lgtseek/lib/";      ### May need to change this depending on where the script is being run
use LGTSeek;
use Bio::Graphics;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use File::Basename;
use run_cmd;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
		'input=s',
		'region=s',
		'sort_bam=s',
		'scale=s',
		'region=s',
		'image_length=s',
		'image_width=s',
		'pad_scale=s',
		'output_prefix=s',
		'output_dir=s',
		'help',
		'help_full'
		);
if($options{help}){die 
"	----------------------------------------------------------------------------------------
	Help: This script will take a bam and use Bio::Graphics to create a .png of the reads.
		----------------------------------------------------------------------------------------
		--input=			Position Sorted bam. Should be very small b/c memory intense. 
		Default output =>		STDOUT. Display - || > foo.png
		--help_full			Complete help information.
		----------------------------------------------------------------------------------------
";
}
if($options{help_full}){ die	
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

my $lgtseek = LGTSeek->new2(\%options);
#############################################################
## Setup & check input
if(!$options{input} && !$ARGV[0]){die "Must give an input: --input=<BAM> or ARGV[0]\n";}
my $input = defined $options{input} ? $options{input} : $ARGV[0];
if($lgtseek->empty_chk({input => $input})==1){die "Error: The input: $options{input} is empty.\n";}
my ($fn,$path,$suff)=fileparse($input,('.srt.bam','.bam'));
#############################################################
## Setup Output
my $output_dir = defined $options{output_dir} ? "$options{output_dir}" : "./";
my $output_prefix =  defined $options{output_prefix} ? "$options{output_prefix}" : "$fn";
my $out = "$output_dir/$output_prefix";
#############################################################
## Setup Defaults
my $sort_bam = defined $options{sort_bam} ? "$options{sort_bam}" : "0";
my $region = defined $options{region} ? "$options{region}" : undef;							## PROBLEM
my $image_width = defined $options{image_width} ? "$options{image_width}" : "1000";
my $image_length = defined $options{image_length} ? "$options{image_length}" : "500";
my $pad_scale = defined $options{pad_scale} ? "$options{pad_scale}" : "0";
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
my $scale_size = defined $options{scale} ? $options{scale} : (($last_position - $first_position)+($pad_scale + $pad_scale));		## = Total + padding on each side
#############################################################
## Setup the graphics panels and tracks.
## Create the Panel
my $panel = Bio::Graphics::Panel->new(
					-length    => $image_length,
					-width     => $image_width,
					-pad_left  => 200,
					-pad_right => 10,
					);
## Scale at the top of the imag
my $scale = Bio::SeqFeature::Generic->new(
	-start => 1,								## Need to make this the start position
	-end   => $scale_size,					## Need to make this the end position
	);
$panel->add_track($scale,
	-glyph   => 'arrow',
	-tick    => 2,
	-fgcolor => 'black',
	-double  => 1,
	);
## Track for bam data 
my $track = $panel->add_track(
					-glyph     => 'graded_segments',
					-label     => 1,
					-bgcolor   => 'blue',
					);

#############################################################
## Build the command to open the bam for reading properly.
my $open_cmd = defined $options{region} ? "samtools view $input $options{region} |" : "samtools view $input |";

#############################################################
## Read through bam and add reads to track
open(BAM,"$open_cmd") || die "Can't open input: $input because: $!\n";
while (<BAM>) {
	chomp;
	my ($id,$flag,$cigar,$start,$sequence) = (split /\t/, $_)[0,1,5,7,9];
	my $converted_flag = $lgtseek->_parseFlag($flag);
	my $length = length $sequence;
	my $end = $start + $length;			## Have to check this!!!
	my $feature = Bio::SeqFeature::Generic->new(
#					-display_name => $id,
					-start        => (($start - $first_position) + $pad_scale),						## +$pad_scale to adjust for padding. shift everything "right" by pad#
					-end          => (($end - $first_position) + $pad_scale),
					);
	$track->add_feature($feature);
}
#############################################################
## Clean up
close BAM;
if($sort_bam==1){run_cmd("rm $output_dir/$fn\_bam2png_tmp_srt.bam");}
if($sort_bam==1){run_cmd("rm $output_dir/$fn\_bam2png_tmp_srt.bai");}
#############################################################

print $panel->png;


__END__


