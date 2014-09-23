#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
use warnings;
use strict;
use Data::Dumper;
use Carp;
use Bio::DB::Fasta;
use Bio::DB::Sam;
use Bio::Graphics;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use GD;
use GD::SVG;
use File::Basename;
use run_cmd;
use print_call;
use bwa;
use LGTSeek;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,           'input|i=s',      'region=s', 'sort|s=i',       'scale=i',       'image_length=i', 'image_width=i',   'pad_scale=i',
    'output_prefix|p=s', 'output_dir|o=s', 'stdout=i', 'draw_stdev|d=i', 'picard_file=s', 'stdev|D=i',      'insert_size|I=i', 'help|?',
) or die "Error: Unrecognized command line option. Please try again.\n";
if ( $options{help} ) { &help; }

my $lgtseek = LGTSeek->new2( \%options );
#############################################################
## Setup & check input
if ( !$options{input} && !$ARGV[0] ) { die "Must give an input: --input=<BAM> or ARGV[0]\n"; }
my $input = defined $options{input} ? $options{input} : $ARGV[0];
my ( $fn, $path, $suff ) = fileparse( $input, ( '.srt.bam', '.bam' ) );
if ( $lgtseek->empty_chk( { input => $input } ) == 1 ) { die "Error: The input: $options{input} is empty.\n"; }
#############################################################
## Setup Output
my $output_dir    = defined $options{output_dir}    ? "$options{output_dir}"    : "$path";
my $output_prefix = defined $options{output_prefix} ? "$options{output_prefix}" : "$fn";
my $out           = "$output_dir/$output_prefix\.svg";
#############################################################
## Setup Defaults
my $sort         = defined $options{sort}         ? $options{sort}         : "0";
my $region       = defined $options{region}       ? $options{region}       : undef;    ## PROBLEM
my $image_width  = defined $options{image_width}  ? $options{image_width}  : "1000";
my $image_length = defined $options{image_length} ? $options{image_length} : "500";
my $pad_scale    = defined $options{pad_scale}    ? $options{pad_scale}    : "0";
#############################################################
## Sort bam
if ( $sort == 1 ) {
    run_cmd("samtools sort $input $output_dir$fn\_bam2svg_tmp_srt");
    run_cmd("samtools index $output_dir$fn\_bam2svg_tmp_srt.bam $output_dir$fn\_bam2svg_tmp_srt.bai");
    $input = "$output_dir$fn\_bam2svg_tmp_srt.bam";
}
#############################################################
## Calculate scale size
my $first_position
    = defined $options{region} ? run_cmd("samtools view -hu $input $options{region} | samtools mpileup -A - | head -n 1 | cut -f2") : run_cmd("samtools mpileup -A $input | head -n 1 | cut -f2");
my $last_position
    = defined $options{region} ? run_cmd("samtools view -hu $input $options{region} | samtools mpileup -A - | tail -n 1 | cut -f2") : run_cmd("samtools mpileup -A $input | tail -n 1 | cut -f2");
my $scale_size = defined $options{scale} ? $options{scale} : ( ( $first_position - $pad_scale ) + ( $last_position + $pad_scale ) );    ## = Total + padding on each side
#############################################################
# Calculate & setup drawing StDev.
my $draw_stdev = defined $options{draw_stdev} ? $options{draw_stdev} : "0";
my $insert_size;
my $stdev;
if ( $draw_stdev == 1 ) {
    if ( defined $options{insert_size} && defined $options{stdev} ) {
        $insert_size = $options{insert_size};
        $stdev       = $options{stdev};
    }
    elsif ( defined $options{picard_file} ) {
        my @lines = `head -n 8 $options{picard_file}`;
        if ( $lines[6] !~ /^MEDIAN_INSERT_SIZE/ ) { confess "Error: The Picard file doesn't look right. Please fix it and try agian.\n"; }
        ( $insert_size, $stdev ) = ( split /\t/, $lines[7] )[ 0, 1 ];
    }
    else { confess "Error: Must use either (--insert_size=[#] & --stdev=[#]) or --picard_file=<file path with Picard insert metrics>\n"; }
}
#############################################################
## Setup the graphics panels and tracks.
## Create the Panel
my $panel = Bio::Graphics::Panel->new(
    -length      => $image_length,
    -width       => $image_width,
    -pad_left    => 10,
    -pad_right   => 10,
    -spacing     => 1,
    -image_class => 'GD::SVG'
);
## Scale at the top of the img
my $scale = Bio::SeqFeature::Generic->new(
    -start => 1,              ## Need to make this the start position
    -end   => $scale_size,    ## Need to make this the end position
);
$panel->add_track(
    $scale,
    -glyph   => 'arrow',
    -tick    => 2,
    -fgcolor => 'black',
    -double  => 1,
);
#############################################################
## Open input Bam data and start to plot it.
my $bam = Bio::DB::Sam->new(
    -bam          => $input,
    -fasta        => $options{ref},
    -expand_flags => 1,
    -autoindex    => 1
);

## Pull read pairs, by region if needed
my @pairs;
if ( $options{region} ) {
    $options{region} =~ m/$(.+)\:(\d+)\-(\d+)^/;
    my ( $seqID, $start, $end ) = ( $1, $2, $3 );
    @pairs = sort( $bam->get_features_by_location(
            -type   => 'read_pair',
            -seq_id => $seqID,
            -start  => $start,
            -end    => $end,
    ) );
}
else {
    @pairs = sort( $bam->features( -type => 'read_pair' ) );
}

## Foreach pair, draw it, and calculate the variance from STDEV + plot
for my $pair (@pairs) {
    my ( $first_mate, $second_mate ) = $pair->get_SeqFeatures;
    my $track;
    if ( $draw_stdev == 1 ) {
        my $read_insert_size = $pair->length;
        my $variance         = $insert_size - $read_insert_size;

        # print STDERR "insert_size: $insert_size\tRead_i-size:".$pair->length."Variance= $variance\n";
        my $color;
        if ( abs($variance) == 0 ) {
            $color = 'white';
        }
        elsif ( abs($variance) <= ( .5 * $stdev ) ) {
            $color = '255,175,175' if $variance < 0;			## Red
            $color = '200,255,200' if $variance > 0;			## Green
        }
        elsif ( abs($variance) <= $stdev ) {
            $color = '255,20,20' if $variance < 0;			## Red
            $color = '115,255,115' if $variance > 0;			## Green
        }
        elsif ( abs($variance) <= ( 2 * $stdev ) ) {
            $color = '185,0,0' if $variance <= 0;			## Red
            $color = '0,165,0' if $variance > 0;			## Green
        }
        elsif ( abs($variance) > ( 2 * $stdev ) ) {
            $color = '100,0,0' if $variance <= 0;			## Red
            $color = '0,75,0' if $variance > 0;			## Green
        }

        # print STDERR "Color: $color\n";
        $track = $panel->add_track(
            -glyph     => 'transcript2',
            -label     => 0,
            -connector => 'dashed',
            -bgcolor   => $color,
        );
    }
    else {
        $track = $panel->add_track(
            -glyph     => 'transcript2',
            -label     => 0,
            -connector => 'dashed',
            -bgcolor   => 'blue',
        );
    }
    $track->add_feature( [ $first_mate, $second_mate ] );
}

#############################################################
my $OFH;
if ( $options{output_dir} || $options{output_prefix} ) {
    open( $OFH, ">", "$out" ) or die "Can't output the output: $out because: $!\n";
}
elsif ( $options{stdout} ) {
    $OFH = *STDOUT;
}
else {
    print STDERR "*** Completed creating the bam.svg, opening display to show results now . . .\n";
    open( $OFH, " | display - " ) or die "Can't output to display.\n";
}
print $OFH $panel->svg;
close $OFH;

if ( $options{sort} ) {
    run_cmd("rm $output_dir$fn\_bam2svg_tmp_srt.bam");
    run_cmd("rm $output_dir$fn\_bam2svg_tmp_srt.bai");
}

sub help {
    die "	----------------------------------------------------------------------------------------
	Help: This script will take a bam and use Bio::Graphics to create a .svg of the reads.
		The bam either needs to cover a small region or --region NEEDS be used. 
		ex: bam2svg.pl foo.bam | display - 
		----------------------------------------------------------------------------------------
		## ARGV[0]			Position sorted bam. Should be very small b/c memory intense. 
		--input=			Position Sorted bam. Should be very small b/c memory intense. 
		  --sort=			<0|1> [0] 1= Sort the input bam on position.
		----------------------------------------------------------------------------------------
		--draw_stdev|d=			<0|1> [0] 1= Color code the reads to illustrate variance from the STDEV of insert size.
						+/- STDEV * .5=Green ; 1=Dark Blue/Red ; 2=Light Blue/Orange
		  --picard_file|P=		< /path/to/file.txt > Picard insert metrics file. Must be used if --stdev && --insert_size are not used.
		  --insert_size|I=		< # > 
		  --stdev|D=			< # >							
		----------------------------------------------------------------------------------------
		--region=			Region of the sorted bam to draw. ex: \"chr19:100-200\" (Must follow this format)
		  --scale=			Size of the scale to draw on the top of the image. [region size]
		  --image_length=		Length of the image. [1000]
		  --image_width= 		Width of the image. [500]
		  --pad_scale=			Number of bp to pad eachside with. [0] ## SHOULD MAKE IT GO NEGATIVE FROM START 
		----------------------------------------------------------------------------------------
		Default output =>	 	Piped into Display
		--stdout=			<0|1> [0] 1= 
		--output_prefix=		[filename] Prefix for the output file. Optional.
		--output_dir=			[input_dir] Directory to put the output. Optional.
		----------------------------------------------------------------------------------------
		--help|?
";
}

__END__

# $track->add_feature(transcript2 => \@pairs); ### doesn't work

### TEST vvvv
	# my $length = $pair->length; 
	# my $f = Bio::Graphics::Feature->new(-segment=>[[$first_mate->start,$first_mate->end],[$second_mate->start,$second_mate->end]],
 #                                  -RGB => 'green');
	# $track->add_feature($f);
	### TEST ^^^^
	# $track->add_feature([$first_mate,$second_mate], -RGB => 'green');

if($pair->length < 300){
		$pair->add_tag_value("RGB","blue");
		# $pair->{RGB}='blue';
		$pair->{attribute}->{RGB}='blue';
		$first_mate->{attribute}->{RGB}='blue';
		$second_mate->{attribute}->{RGB}='blue';
		# $pair->{RGB}=['blue'];
		# $pair->{tag}->{RGB}=['blue'];
		# $first_mate->{tag}->{RGB}=['blue'];
		# $second_mate->{tag}->{RGB}=['blue'];
		# $first_mate->{RGB}=['blue'];
		# $second_mate->{RGB}=['blue'];
	} else {
		$pair->add_tag_value("RGB","red");
		# $pair->{RGB}='red';
		# $first_mate->{RGB}='red';
		# $pair->{RGB}=['red'];
		# $first_mate->{RGB}=['red'];
		# $second_mate->{RGB}='red';
		$pair->{attribute}->{RGB}='red';
		$first_mate->{attribute}->{RGB}='red';
		$second_mate->{attribute}->{RGB}='red';
		# $first_mate->{tag}->{RGB}=['red'];
		# $second_mate->{tag}->{RGB}=['red'];
	}
	print STDERR "*** DUMPER ***\n";
	print Dumper($pair);
	

	#### TEST vvvv
	my $insert_size = 256;
	my $stdev = 20;
	my $read_i_size = $length;
	print STDERR "Read-insert-size: $read_i_size\n";
	my $variance = $insert_size - $read_i_size;
	print STDERR "$variance\n";
	if (abs($variance) < (.5 * $stdev)){
		$color = 'green';
	} elsif (abs($variance) <= $stdev){
		$color = 'red' if $variance < 0;
		$color = 'blue' if $variance > 0;
	} elsif (abs($variance) <= (2 * $stdev)){
		$color = 'orange' if $variance < 0;
		$color = 'cyan' if $variance > 0;
	}
	print STDERR "$color\n";
	#### TEST ^^^^

	$first_mate->{sam}->{RGB}=[$color];
	$second_mate->{sam}->{RGB}=[$color];
	print Dumper($first_mate);
	
	my $length = $pair->length; 

	#-bgcolor => "featureRGB",
	#-RGB => $color,

sub {
									my $feature = shift;
									# print Dumper($feature);
									# my @mates = $feature->(-type => 'read_pair');
									print STDERR "MATES:\n";
									print Dumper(@mates);

									my $insert_size = 156;
									my $stdev = 20;
									my $read_i_size = $feature->length;
									print STDERR "Read-insert-size: $read_i_size\n";
									my $variance = $insert_size - $read_i_size;
									if (abs($variance) < (.5 * $stdev)){
										return 'green';
									} elsif (abs($variance) <= $stdev){
										return 'red' if $variance < 0;
										return 'blue' if $variance > 0;
									} elsif (abs($variance) <= (2 * $stdev)){
										return 'orange' if $variance < 0;
										return 'cyan' if $variance > 0;
									}
								}
