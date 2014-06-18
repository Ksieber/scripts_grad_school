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
my $results = GetOptions (\%options,
	'bam1=s',
	'bam2=s',
	'ref1=s',
	'ref2=s',
	'bam1_region=s',
	'bam2_region=s',
	'ref1_region=s',
	'ref2_region=s',
	'sort1=i',
	'sort2=i',
	'reads_list=s',
	'n_num|n=i',
	'draw_nstring=i',
	'fix_orientation=i',
	'M_only=i',
	'MM_only=i',
	'png=i',
	'svg=i',
	'draw_stdev|d=i',
	'draw_both|B=i',
	'picard_file|P=s',
	'stdev|D=i',
	'insert_size|I=i',
	'image_length=i',
	'image_width=i',
	'pad_scale=i',
	'output_dir|o=s',
	'output_prefix|p=s',
	'merged_ref_name=s',
	'help|h'
) or confess "Invalid arguement. Please try agian.\n";

if($options{help}){&help;}   ## &help is @ the end of the script
if(!$options{ref1} || !$options{ref2} || !$options{bam1}) {
	confess "Error: You MUST pass ATLEAST the folowing arguements: --bam1 --ref1 --ref2 . Please try agian.\n";
}

print_call(\%options);
## Setup a few Default values
my $lgtseek = LGTSeek->new2(\%options);
my $fix_orientation = defined $options{fix_orientation} ? $options{fix_orientation} : "1";
my $draw_png = defined $options{draw_png} ? $options{draw_png} : "0";

## Get filenames, paths, create output directory and set output names.
my ($fn1,$path1,$suff1) = fileparse($options{bam1},qr/\.[^\.]+/);
my $out_dir = $options{output_dir} ? $options{output_dir} : $path1;
run_cmd("mkdir -p $out_dir");
my $output_prefix = $options{output_prefix} ? $options{output_prefix} : $fn1;
if($output_prefix=~/(.+)\_psort$/){
	$output_prefix = $1;		## Remove a trailing _psort form input prefix 
}
my $out = $out_dir . $output_prefix;
my $working_dir = "$out_dir"."/merge_2lgt_bams/";
run_cmd("mkdir -p $working_dir"); 		# Create a safe directory to work in
my $log = "$working_dir\/log.txt";
run_cmd("touch $log");
print_notebook(\%options);

## Map bam1 against ref2 if bam2 wasn't given.
if(!$options{bam2} && $options{ref2}){
	my ($fnR,$pathR,$suffR) = fileparse($options{ref2},qr/\.[^\.]+/);
	$options{bam2} = &bwa_align($options{bam1},$options{ref2},{output_prefix => "$fn1\_$fnR", output_dir => $working_dir, sort_index => 1, cmd_log => 1});
}

my $sort1 = defined $options{sort1} ? $options{sort1} : 0;
my $sort2 = defined $options{sort2} ? $options{sort2} : 0;

## Comment in our out according to samtools version
## Use this with samtools version >= 0.2.0 
if($sort1==1){
	run_cmd("samtools sort -o $options{bam1} - | samtools view -b - > $working_dir/$fn1\_psort.bam",$log);
	run_cmd("samtools index $working_dir/$fn1\_psort.bam",$log);
	$options{bam1} = "$working_dir/$fn1\_psort.bam";
}
if($sort2==1){
	my ($fn2,$path2,$suff2) = fileparse($options{bam2},qr/\.[^\.]+/); 
	run_cmd("samtools sort -o $options{bam2} - | samtools view -b - > $working_dir\/$fn2\_psort.bam",$log);
	run_cmd("samtools index $working_dir\/$fn2\_psort.bam",$log);
	$options{bam2} = "$working_dir\/$fn2\_psort.bam";
}

my %reads;  ## Hash of desired read id's.
if($options{reads_list}){
	open(IN,"<","$options{reads_list}") or confess "Can't open reads_list: $options{reads_list} because: $!\n";
	while(<IN>){
		chomp; 
		$reads{$_}++;
	}
	close IN;
}

### Process files
# Grab important bam data
my $bam1_data = &pull_bam_data({ 
	bam 	=> $options{bam1},
	region 	=> $options{bam1_region} 
	});
my $bam2_data = undef;
if($options{bam2}){
	$bam2_data = &pull_bam_data({ 
		bam 	=> $options{bam2},
		region 	=> $options{bam2_region} 
		});
}

my $merged_ids = &merge_hash_ids( [$bam1_data->{id_hash}, $bam2_data->{id_hash}] );

# Grab ref seq. for each
my $ref1_region = &pull_ref_data({
	ref => $options{ref1}, 
	ref_region => $options{ref1_region},
	bam_data => $bam1_data, 
	merged_ids => $merged_ids,	
	});

my $ref2_region = undef;
if($options{bam2}){
	$ref2_region = &pull_ref_data({
		ref => $options{ref2}, 
		ref_region => $options{ref2_region},
		bam_data => $bam2_data, 
		merged_ids => $merged_ids,
		});
}

my $merged_ref = &merge_refs({
	bam1_data 	=> 	$bam1_data,
	bam2_data 	=> 	$bam2_data,
	ref1_region =>	$ref1_region,
	ref2_region	=>	$ref2_region,
	output_dir	=> $out_dir,
	});

## Merge the sam reads into a bam only for mapping @ the new reference.
# Create a fake header for the tmp. mapping bam.
open(SAM,">","$working_dir/Merged-reads-only.sam") or confess "Error: Couldn't open output merged.sam: $working_dir/Merged-reads-only.sam";
print SAM @{$bam1_data->{header}};
print SAM @{$bam2_data->{header}} if($options{bam2});

# Pring the reads to the tmp mappig bam by using id's
foreach my $id (keys %{$merged_ids}){
	print SAM "$bam1_data->{id_hash}->{$id}\n";
	print SAM "$bam2_data->{id_hash}->{$id}\n" if ($options{bam2});
}
close SAM;

# Convert the tmp mapping sam -> bam.
run_cmd("samtools view -S $working_dir/Merged-reads-only.sam -bo $working_dir\/Merged-reads-only.bam",$log);
run_cmd("rm $working_dir\/Merged-reads-only.sam",$log);

# Map the tmp bam @ the new reference.
if($options{M_only}){
	&bwa_align("$working_dir\/Merged-reads-only.bam",$merged_ref->{file},{output_prefix => "Merged\-final\_$output_prefix", output_dir => $out_dir, M_only => $options{M_only}, cmd_log => 1}); 
} elsif ($options{MM_only}) { 
	&bwa_align("$working_dir\/Merged-reads-only.bam",$merged_ref->{file},{output_prefix => "Merged\-final\_$output_prefix", output_dir => $out_dir, MM_only => $options{MM_only}, cmd_log => 1});
} else {
	&bwa_align("$working_dir\/Merged-reads-only.bam",$merged_ref->{file},{output_prefix => "Merged\-final\_$output_prefix", output_dir => $out_dir, cmd_log => 1});
}
run_cmd("cat $out_dir/log.txt >> $log",$log);
run_cmd("rm $out_dir/log.txt",$log);
# Make a position sorted final output bam in the /working/dir/ to be used to draw the img
run_cmd("samtools sort $out_dir\/Merged\-final\_$output_prefix\.bam $working_dir\/Merged\-final\_$output_prefix\-psrt",$log);
my $psrt_bam = "$working_dir\/Merged\-final\_$output_prefix\-psrt.bam";

## Iterate to draw both STDEV & Normal img's if --draw_both
if($options{draw_both}){
	## First draw the "normal" img
	&draw_img({ 
		bam => $psrt_bam, 
		ref_data => $merged_ref, 
		output_dir => $out_dir,
		output_prefix => "Merged\-$output_prefix",
		png => $options{png}, 
		svg => $options{svg},
		draw_stdev => "0",
		draw_nstring => $options{draw_nstring},
		});
	## Second draw the color coded img
	&draw_img({ 
		bam => $psrt_bam, 
		ref_data => $merged_ref,
		output_dir => $out_dir,
		output_prefix => "Merged\-$output_prefix\_stdev", 
		png => $options{png}, 
		svg => $options{svg},
		draw_stdev => "1",
		stdev => $options{stdev},
		insert_size => $options{insert_size},
		picard_file => $options{picard_file},
		draw_nstring => $options{draw_nstring},
		});
} elsif($options{png} || $options{svg}){ 
	&draw_img({ 
		bam => $psrt_bam, 
		ref_data => $merged_ref,
		output_dir => $out_dir,
		output_prefix => "Merged\-final\-$output_prefix", 
		png => $options{png}, 
		svg => $options{svg},
		draw_stdev => $options{draw_stdev},
		stdev => $options{stdev},
		insert_size => $options{insert_size},
		picard_file => $options{picard_file},
		draw_nstring => $options{draw_nstring},
		});
}

print_complete(\%options);
## Done processing

###################################
########### Subroutines ###########
###################################

=head1

Title   : pull_bam_data
Usage   : my $bam_data = pull_bam_data($bam,$region);
Function: Extract the reads from a specific bam region and return them in an array. 
Args 	: 
		bam 		=> /file/path/input.bam
		region 		=> 'chr:100-200'
Returns : Returns a data structure with an array of the data as well as the data stored by hash.
   	 	'hash'			=> \%bam_data_hash,
   	 	'header'		=> \@samtools_header
		'strand'		=> Integer. Positive means more reads map to + strand.
=cut
sub pull_bam_data {
	my $opts = shift;
	if($lgtseek->empty_chk({input => $opts->{bam}})==1){die "Error: The input: $opts->{bam} is empty.\n";}
	
	# Grab the samtools header, and remove the @PG line.
	my @header = `samtools view -H $opts->{bam} 2>>$log` or confess "Error: Couldn't get the samtools header from: $opts->{bam}\n";
	delete $header[-1];

	my $open_cmd = defined $opts->{region} ? "samtools view -F0x4 $opts->{bam} \"$opts->{region}\" 2>>$log |" : "samtools view -F0x4 $opts->{bam} 2>>$log |";

	my $bam_orientation = 0;
	my %bam_data;
	open(BAM,"$open_cmd") or confess "Error: Unable to open input bam: $opts->{bam} because: $!\n";
	while(<BAM>){
		chomp;
		my @f=split;
		## If a read list was passed, make sure this is a read we want
		my $read_id = $f[0];
		if($options{reads_list}){ next if(!$reads{$read_id}); }
		## Add bam data to hash by id
		$bam_data{$read_id}=$_;
	}
	close BAM;

    my $retval = {
        'id_hash' => \%bam_data,
        'header' => \@header,
        'strand' => "0",
    };
	return $retval;
}


=head1

Title   : pull_ref_data
Usage   : my $ref_seq = pull_ref_data($bam_data_in_array,$ref,$region);
Function: Extract the reference sequence from a specific fasta region
Args 	: 
		ref 		=>	{Mandatory}{Priority}
		ref_region 	=>	{Optional}{Priority}
		bam_data 	=> 	{Mandatory}
		merged_ids 	=>	{Mandatory}
Returns : Bio::DB::Fasta->seq obj. 

=cut

sub pull_ref_data {
	my $opts = shift;

	if($lgtseek->empty_chk({input => $opts->{ref}})==1){die "Error: The ref: $opts->{ref} is empty.\n";}

	my $ref_lower_range;
	my $ref_upper_range;
	my $ref_chr;
	my $bam_lower_range;
	my $bam_upper_range;
	my $bam_chr;

	if($opts->{ref_region}){
		$opts->{ref_region} =~ /^([\w\-\|\.]+)\:(\d+)\-(\d+)$/;
		($ref_chr,$ref_lower_range,$ref_upper_range)= ($1,$2,$3);
	} else {
		my $first_read = 0;
		foreach my $read_id (keys %{$opts->{merged_ids}}){
			if($first_read==0){
				my $bam_line = $opts->{bam_data}->{id_hash}->{$read_id};
				my ($chr,$position,$seq) = (split/\t/,$bam_line)[2,3,9];
				$bam_chr = $chr;
				$bam_lower_range = $position;
				$bam_upper_range = $position + length($seq);
				$first_read++;
			} else {
				my $bam_line = $opts->{bam_data}->{id_hash}->{$read_id};
				my ($chr,$position,$seq) = (split/\t/,$bam_line)[2,3,9];
				next if($chr ne $bam_chr);
				$bam_lower_range = $position if ($position < $bam_lower_range);
				$bam_upper_range = ($position + length($seq)) if (($position + length($seq)) > $bam_upper_range);
			}
		}
	}

	# Use reference region is it exist. 
	my $actual_lower_range = defined $opts->{ref_region} ? $ref_lower_range : $bam_lower_range;
	my $actual_upper_range = defined $opts->{ref_region} ? $ref_upper_range : $bam_upper_range;	
	my $actual_chr = defined $opts->{ref_region} ? $ref_chr : $bam_chr; 

	# Finally, pull the region from the fasta reference
	my $ref_db = Bio::DB::Fasta->new($opts->{ref});
	my $seq = $ref_db->seq($actual_chr,$actual_lower_range => $actual_upper_range);
	
	my $retval = {
		'seq' => $seq,
		'range' => "$actual_chr\:$actual_lower_range\-$actual_upper_range",
	};
	return $retval;
}

=head1

Title   : create_n_string
Usage   : my $n_string = &create_n_string($n_num);
Function: Create a string of N's to insert inbetween merged references. 
Returns : String of N's

=cut
sub create_n_string {
	my $n_num = shift;
	my $n_string;
	for (my $i=0; $i < $n_num; $i++){
		$n_string = $n_string . "N";
	}
	return $n_string;
}

=head1

Title   : merge_hash_ids
Usage   : my $merged_ids = create_n_string([%hash1,%hash2,%hash3]);
Function: Merge keys from multiple hashes to find keys that exist in all hashes
Returns : 1 hash with uniq id's

=cut
sub merge_hash_ids {
	my $ref_array_of_hashes = shift;
	my %seen_ids;
	my $number_of_hashes = @{$ref_array_of_hashes};
	foreach my $hash (@{$ref_array_of_hashes}){
		foreach my $ids (keys %{$hash}){
			$seen_ids{$ids}++;
		}
	}
	my %merged_ids;
	foreach my $id (keys %seen_ids){
		if($seen_ids{$id}==$number_of_hashes){
			$merged_ids{$id}++;
		}
	}
	return \%merged_ids;
}

=head1

Title   : merge_refs
Usage   : my 
Function: 
Args    : 
		bam1_data 	=> 	$bam1_data
		bam2_data 	=> 	$bam2_data
		ref1_region =>	$ref1_region
		ref2_region	=>	$ref2_region
		output_dir	=>  $output_dir (for ref_ranges.txt)

Returns : 1 hash with uniq id's

=cut
sub merge_refs {
	my $opts = shift; 

	# Determine & create the string of N's to put between the references for merging. 
	my $n_num = defined $options{n_num} ? $options{n_num} : 100; 
	my $n_string = &create_n_string($n_num);
	
	my $ref1_region;
	my $ref2_region;

	## Fix orientation based on majority of the sequencing reads if the options is passed
	if($fix_orientation==1){ 
		map {
			my $bam_data = $_;
			foreach my $id (keys %{$merged_ids}){
				my $bam_line = $bam_data->{id_hash}->{$id};
				my $raw_flag = (split/\t/,$bam_line)[1];
				my $parsed_flag = $lgtseek->_parseFlag($raw_flag);
				if($parsed_flag->{'qrev'}){
					$bam_data->{strand}--;
				} else {
					$bam_data->{strand}++;
				}
			}
		} ($opts->{bam1_data},$opts->{bam2_data});
		if($opts->{bam1_data}->{strand}<=0 && $opts->{bam2_data}->{strand}>=0){
			$ref1_region = $opts->{ref2_region};
			$ref2_region = $opts->{ref1_region};
		}
	} else {
		$ref1_region = $opts->{ref1_region};
		$ref2_region = $opts->{ref2_region};
	}


	## Build the new reference Sequence
	my $new_ref = $ref1_region->{seq} . $n_string . $ref2_region->{seq};
	
	## Calculate the cordinates for the new reference
	my $cords = {
		ref1_lower 	=> 1,
		ref1_upper 	=> length($ref1_region->{seq}),
		n_lower 	=> length($ref1_region->{seq})+1,
		n_upper 	=> length($ref1_region->{seq})+length($n_string),
		ref2_lower 	=> length($ref1_region->{seq})+length($n_string)+1,
		ref2_upper 	=> length($ref1_region->{seq})+length($n_string)+length($ref2_region->{seq}),
	};

	## Open text file to print the cordinates of the references and the cordinates of the new image. 
	open (TXT,">","$opts->{output_dir}/Ref_img_cords.txt") or confess "Error: can't open the output file for the ref & img cordinates: $opts->{output_dir}/Ref_img_cords.txt because: $!";
	print TXT "Reference\tImage_cords\n";
	print TXT "$ref1_region->{range}\t$cords->{ref1_lower}\-$cords->{ref1_upper}\n";
	print TXT "n_string\:1-".length($n_string)."\t$cords->{n_lower}\-$cords->{n_upper}\n";
	print TXT "$ref2_region->{range}\t$cords->{ref2_lower}\-$cords->{ref2_upper}\n";
	close TXT; 

	## Print new reference Sequence.
	my $new_ref_name = defined $options{merged_ref_name} ? $options{merged_ref_name} : "Merged";
	open (FASTA,">","$working_dir\/Merged-refs\_$output_prefix.fa") or confess "Error: can't open the output file for the new merged fasta: $working_dir\/Merged-tmp\_$output_prefix.fa because: $!";
	print FASTA "\>$new_ref_name\n";	## Print new "chr" header line for .fasta
	print FASTA $new_ref . "\n";
	close FASTA;
	run_cmd("bwa index $working_dir\/Merged-refs\_$output_prefix.fa 2>>$log",$log);

	return {
		file 		=> "$working_dir\/Merged-refs\_$output_prefix.fa",
		seq 		=> $new_ref,
		cords		=> $cords,
	};

}

sub draw_img {
	my $opts = shift;
	my $bam = $opts->{bam};
	my $merged_ref = $opts->{ref_data};
	my $png = defined $opts->{png} ? $opts->{png} : "0";
	my $svg = defined $opts->{svg} ? $opts->{svg} : "0";
	my $draw_stdev = defined $opts->{draw_stdev} ? $opts->{draw_stdev} : "0";
	
	if($lgtseek->empty_chk({input => $bam})==1){confess "Error: The Merged-final.bam: $bam is empty.\n";}
	#############################################################
	## Setup OUTPUT
	my ($fn,$path,$suff)=fileparse($bam,'.bam');
	my $out_dir = defined $opts->{output_dir} ? $opts->{output_dir} : $path;
	my $out_pref = defined $opts->{output_prefix} ? $opts->{output_prefix} : $fn;
	my $out = "$out_dir\/$out_pref";
	#############################################################
	## Setupt image boundries
	my $image_width = defined $options{image_width} ? $options{image_width} : "1600";
	my $image_length = defined $options{image_length} ? $options{image_length} : "800";
	my $pad_scale = defined $options{pad_scale} ? $options{pad_scale} : "0";
	#############################################################
	## Calculate scale size ##
	my $scale_size = defined $options{scale} ? $options{scale} : length($merged_ref->{seq});
	#############################################################
	# Calculate & setup drawing StDev.
	my $insert_size;
	my $stdev; 
	if($draw_stdev==1){
		if(defined $opts->{insert_size} && defined $opts->{stdev}){
			$insert_size = $opts->{insert_size};
			$stdev = $opts->{stdev};
		} elsif (defined $opts->{picard_file}){
			my @lines = `head -n 8 $opts->{picard_file}`;
			if($lines[6]!~/^MEDIAN_INSERT_SIZE/){ confess "Error: The Picard file doesn't look right. Please fix it and try agian.\n";}
			($insert_size,$stdev) = (split/\t/,$lines[7])[0,1];
		} else {
			confess "Error: Must use either (--insert_size=[#] & --stdev=[#]) or --picard_file=<file path with Picard insert metrics>\n";
		}
	}
	#############################################################
	# Initialize Bio:Graphics img
	my @panel_options = (
		-length    => $image_length,
		-width     => $image_width,
		-pad_left  => 10,
		-pad_right => 10,
		-spacing   => 1,
		);
	if($svg==1){push(@panel_options,(-image_class=>'GD::SVG'))};
	my $panel = Bio::Graphics::Panel->new(@panel_options);
	#############################################################
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
	#############################################################
	## Tracks for references
	## Left Reference
	my $ref1_feature = Bio::SeqFeature::Generic->new(
		-start 		=> $merged_ref->{cords}->{ref1_lower},
		-end 		=> $merged_ref->{cords}->{ref1_upper},
		);
	$panel->add_track($ref1_feature,
		-glyph 		=> 'generic',
		-bgcolor 	=> 'red',
		);
	# ## N-string "contig"
	if($opts->{draw_nstring}){
		my $nstring_feature = Bio::SeqFeature::Generic->new(
			-start 		=> $merged_ref->{cords}->{n_lower},
			-end 		=> $merged_ref->{cords}->{n_upper},
			);
		$panel->add_track($nstring_feature,
			-glyph 		=> 'crossbox',
			-fgcolor 	=> 'black',
			-bgcolor	=> '128,128,128'
			);
	}
	## Right Reference
	my $ref2_feature = Bio::SeqFeature::Generic->new(
		-start 		=> $merged_ref->{cords}->{ref2_lower},
		-end 		=> $merged_ref->{cords}->{ref2_upper},
		);
	$panel->add_track($ref2_feature,
		-glyph 		=> 'generic',
		-bgcolor 	=> 'blue',
		);
	#############################################################
	## Track for bam data 
	## Initialize the bam db
	my $sam = Bio::DB::Sam->new(
		-bam => $bam,
		-fasta => $merged_ref->{file},
		-expand_flags => "1",
		-autoindex => "1",
		);
	## Add feature by reads_list and region later
	my @pairs = sort($sam->features(-type => 'read_pair'));
	for my $pair (@pairs){
		my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
		my $track;
		if($draw_stdev==1){
			my $read_insert_size = $pair->length;
			my $variance = $insert_size - $read_insert_size;
			my $color;
			if (abs($variance)==0){
				$color = 'white';
			} elsif (abs($variance) <= (.5 * $stdev)){
				$color = '255,234,234' if $variance < 0;
				$color = '234,255,234' if $variance > 0;
			} elsif (abs($variance) <= $stdev){
				$color = '255,140,140' if $variance < 0;
				$color = '140,255,140' if $variance > 0;
			} elsif (abs($variance) <= (2 * $stdev)){
				$color = '255,0,0' if $variance <= 0;
				$color = '0,255,0' if $variance > 0;
			} elsif (abs($variance) > (2 * $stdev)){
				$color = '170,0,0' if $variance <= 0;
				$color = '0,170,0' if $variance > 0;
			}
			$track = $panel->add_track(
							-glyph     => 'transcript2',
							-label     => 0,
							-connector => 'dashed',
							-bgcolor   => $color,
							);
		} else {
			$track = $panel->add_track(
							-glyph     => 'transcript2',
							-label     => 0,
							-connector => 'dashed',
							-bgcolor   => 'purple',
							);
		}
		$track->add_feature([$first_mate,$second_mate]);
	}

	#############################################################
	## Print out the image
	my $OFH;
	my $output_suffix = $png==1 ? ".png" : ".svg";
	if($opts->{output_dir} || $opts->{output_prefix}){
		open($OFH,">","$out$output_suffix") or die "Can't output the output: $out-Merged.png because: $!\n";
	} elsif($opts->{stdout}){
		$OFH = *STDOUT;
	} else {
		open ($OFH, " | display - ") or die "Can't output to display.\n";
	}
	print $OFH $panel->png if $png==1;
	print $OFH $panel->svg if $svg==1;
	close $OFH;
}

sub help {
	die "Help: This script will take a 2 bams mapped @ different references and merge the references & map at the merged reference. The final new bam can be drawn as a png or svg. This is useful for merging and illustrating LGT regions.
		------------------------------------------------------------------------------------------------------------------------------------------------------------
		--bam1=				bam1. Assumes position sorted & indexed. If not, use --sort=1. (Mandatory)
		  --bam1_region=		<chr#:100-200> Pull reads only from this region. (Highly recommended)
		  --sort1=			<0|1> [0] 1= Position sort & index bam1.
		--bam2=				bam2. Assumes position sorted & indexed. If not, use --sort=1. (Optional)
		  --bam2_region=		<chr#:100-200> Pull reads only from this region. (Highly recommended)
		  --sort2=			<0|1> [0] 1= Position sort & index bam2.
		------------------------------------------------------------------------------------------------------------------------------------------------------------
		--ref1=				Reference 1 fasta. 		(Assumes bam1 is already mapped aginst ref1)(Mandatory)
		  --ref1_region=		<chr#:100-200> Use this reference range to map & draw reads against. If no ref_region, the bam reads' region will be used. 
		--ref2=				Reference 2 fasta.		(If no --bam2, and --ref2 is used, --bam1 will be mapped against --ref2)
		  --ref2_region=		<chr#:100-200> Use this reference range to map & draw reads against. If no ref_region, the bam reads' region will be used. 
		------------------------------------------------------------------------------------------------------------------------------------------------------------
		--reads_list=			Path to a file with a list of desired reads to parse for. 1 read ID / line. 
		--M_only=			<0|1> [0] 1= When remapping to the merged reference, only keep M_* read pairs
		--MM_only=    			<0|1> [1] 1= When remapping to the merged reference, only keep M_M read pairs (Highly recommended)
		--merged_ref_name= 		Name for the new reference.  		[Merged]
		--n_num|n= 			Number of \"N's\" to insert inbetween merged references. [100]
		--draw_nstring=	<0|1> [0] 1= Draw the n-string \"contig\".
		------------------------------------------------------------------------------------------------------------------------------------------------------------
		--png=				<0|1> [0] 1=Create a png img of the merged bam.
		--svg=				<0|1> [0] 1=Create a svg img of the merged bam.
		  --image_length=		Ajust the length of the png created.
		  --image_width=		Adjust the width of the png created. 
		  --pad_scale=			Pad white space around img. 
		  --fix_orientation=		<0|1> [1] 1= Try to determine how the references should be organized L-vs-R to make Mates face eachother.
							  0= Bam1 is on Left, Bam2 is on Right.
		------------------------------------------------------------------------------------------------------------------------------------------------------------
		--draw_stdev|d=			<0|1> [0] 1=Color code the reads based on variance from STDEV of insert size;
						+/- STDEV * .5=Green ; 1=Dark Blue/Red ; 2=Light Blue/Orange
		--draw_both|B=			<0|1> [0] 1= Draw both a \"normal\" & stdev color coded img.
		  --picard_file|P=		< /path/to/file.txt > Picard insert metrics file. Must be used if --stdev && --insert_size are not used.
		  --insert_size|I=		< # > 
		  --stdev|D=			< # >
		------------------------------------------------------------------------------------------------------------------------------------------------------------							  
		--output_dir|o=			Directory for output. 				[/options/bam1/dir/]
		--output_prefix|p=		Prefix for the output fasta & bam. 	[bam1-merged-bam2]
		--stdout=				<0|1> [0] 1= Output goes to STDOUT. Either pipe it into a \"display\" ( | display - ) or redirect it to a new file ( > new.img)
		--help|h
		------------------------------------------------------------------------------------------------------------------------------------------------------------
		Example: perl merge_2lgt_bams.pl --bam1=bam_with_lgt_reads.bam --sort1=1 --ref1=hg19.fa --ref2=bacteria_with_lgt.fa --ref1_region=chr1:100-500 --bam1_region=chr1:200-400 --MM_only=1 --svg=1
				This will pull reads from the {bam_with_lgt_reads.bam} @ {chr1:200-400} ; map them against {bacteria_with_lgt.fa} ; create a new reference from {chr1:100-500} & 
				the region where the reads mapped in the {bacteria_with_lgt.fa}. The reads that are {mapped-mapped} to the new reference are then kept in the final bam and drawn as an {svg}. 
		------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
}

###################################

