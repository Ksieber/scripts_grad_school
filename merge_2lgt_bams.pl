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
	'fix_orientation=i',
	'M_only=i',
	'MM_only=i',
	'png=i',
	'svg=i',
	'image_length=i',
	'image_width=i',
	'pad_scale=i',
	'output_dir|o=s',
	'output_prefix=s',
	'merged_ref_name=s',
	'help|h'
) or confess "Invalid arguement. Please try agian.\n";

if($options{help}){&help;}   ## &help is @ the end of the script
if(!$options{ref1} || !$options{ref2} || !$options{bam1}) {
	confess "Error: You MUST pass ALL of the folowing arguements: --bam1 --ref1 --ref2 . Please try agian.\n";
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
my $out = $out_dir . $output_prefix;
my $working_dir = "$out_dir"."/merge_2lgt_bams/";
run_cmd("mkdir -p $working_dir"); 		# Create a safe directory to work in
my $log = "$working_dir\/log.txt";
run_cmd("touch $log");

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

## Use this with samtools version <= 0.1.19
# if($sort1==1){
# 	run_cmd("samtools sort $options{bam1} $working_dir/$fn1\_psort",$log);
# 	run_cmd("samtools index $working_dir/$fn1\_psort.bam",$log);
# 	$options{bam1} = "$working_dir/$fn1\_psort.bam";
# }
# if($sort2==1){
# 	my ($fn2,$path2,$suff2) = fileparse($options{bam2},qr/\.[^\.]+/); 
# 	run_cmd("samtools sort $options{bam2} $working_dir\/$fn2\_psort",$log);
# 	run_cmd("samtools index $working_dir\/$fn2\_psort.bam",$log);
# 	$options{bam2} = "$working_dir\/$fn2\_psort.bam";
# }

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

if($options{png} || $options{svg}){ 
	&draw_img({ 
		bam => "$out_dir\/Merged\-final\_$output_prefix\.bam", 
		ref_data => $merged_ref, 
		png => $options{png}, 
		svg => $options{svg} 
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
	return $seq;
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

Title   : merge_hash_ids
Usage   : my $merged_ids = create_n_string([%hash1,%hash2,%hash3]);
Function: Merge keys from multiple hashes to find keys that exist in all hashes
Args    : 
		bam1_data 	=> 	$bam1_data
		bam2_data 	=> 	$bam2_data
		ref1_region =>	$ref1_region
		ref2_region	=>	$ref2_region

Returns : 1 hash with uniq id's

=cut
sub merge_refs {
	my $opts = shift; 

	# Determine & create the string of N's to put between the references for merging. 
	my $n_num = defined $options{n_num} ? $options{n_num} : 100; 
	my $n_string = &create_n_string($n_num);
	
	# Determine how the two references should be facing
	my $new_ref = $opts->{ref1_region} . $n_string . $opts->{ref2_region};
	## Cords = Ref1 : nstring : ref2
	my $cords = {
					ref1_lower 	=> 1,
					ref1_upper 	=> length($opts->{ref1_region})+1,
					n_lower 	=> length($opts->{ref1_region})+2,
					n_upper 	=> length($opts->{ref1_region})+2+length($n_string),
					ref2_lower 	=> length($opts->{ref1_region})+2+length($n_string)+1,
					ref2_upper 	=> length($opts->{ref1_region})+2+length($n_string)+1+length($opts->{ref2_region}),
				};

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
			$new_ref = $opts->{ref2_region} . $n_string . $opts->{ref1_region};
			$cords = {
					ref1_lower 	=> 1,
					ref1_upper 	=> length($opts->{ref2_region})+1,
					n_lower 	=> length($opts->{ref2_region})+2,
					n_upper 	=> length($opts->{ref2_region})+2+length($n_string),
					ref2_lower 	=> length($opts->{ref2_region})+2+length($n_string)+1,
					ref2_upper 	=> length($opts->{ref2_region})+2+length($n_string)+1+length($opts->{ref1_region}),
					};
		}
	} 

	my $new_ref_name = defined $options{merged_ref_name} ? $options{merged_ref_name} : "Merged";
	open(FASTA,">","$working_dir\/Merged-refs\_$output_prefix.fa") or confess "Error: can't open the output file for the new merged fasta: $working_dir\/Merged-tmp\_$output_prefix.fa because: $!";
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

	if($lgtseek->empty_chk({input => $bam})==1){confess "Error: The Merged-final.bam: $bam is empty.\n";}

	my ($fn,$path,$suff)=fileparse($bam,'.bam');
	run_cmd("samtools sort $bam $working_dir/$fn\-psrt",$log);
	my $psrt_bam = "$working_dir/$fn\-psrt.bam";
	my $sam = Bio::DB::Sam->new(
		-bam => $psrt_bam,
		-fasta => $merged_ref->{file},
		-expand_flags => 1,
		-autoindex => 1
		);
	#############################################################
	my $image_width = defined $options{image_width} ? $options{image_width} : "1600";
	my $image_length = defined $options{image_length} ? $options{image_length} : "800";
	my $pad_scale = defined $options{pad_scale} ? $options{pad_scale} : "0";
	#############################################################
	## Calculate scale size ##
	my $scale_size = defined $options{scale} ? $options{scale} : length($merged_ref->{seq});
	#############################################################
	# Initialize Bio:Graphics img
	my @panel_options = (
		-length    => $image_length,
		-width     => $image_width,
		-pad_left  => 10,
		-pad_right => 10,
		);
	if($svg==1){push(@panel_options,(-image_class=>'GD::SVG'))};
	my $panel = Bio::Graphics::Panel->new(@panel_options);

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

	## Track for references
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
	# my $nstring_feature = Bio::SeqFeature::Generic->new(
	# 	-start 		=> $merged_ref->{cords}->{n_lower},
	# 	-end 		=> $merged_ref->{cords}->{n_upper},
	# 	);
	# $panel->add_track($nstring_feature,
	# 	-glyph 		=> 'crossbox',
	# 	-fgcolor 	=> 'black',
	# 	);
	## Right Reference
	my $ref2_feature = Bio::SeqFeature::Generic->new(
		-start 		=> $merged_ref->{cords}->{ref2_lower},
		-end 		=> $merged_ref->{cords}->{ref2_upper},
		);
	$panel->add_track($ref2_feature,
		-glyph 		=> 'generic',
		-bgcolor 	=> 'blue',
		);
		
	## Track for bam data 
	my $bam_track = $panel->add_track(
		-glyph     => 'transcript2',
		-label     => 0,
		-bgcolor   => 'purple',
		);

	## Add feature by reads_list and region later
	my @pairs = $sam->features(-type => 'read_pair');
	for my $pair (@pairs){
		my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
		$bam_track->add_feature([$first_mate,$second_mate]);
	}

	#############################################################
	my $OFH;
	my $output_suffix = $png==1 ? ".png" : ".svg";
	if($options{output_dir} || $options{output_prefix}){
		open($OFH,">","$path\/$fn$output_suffix") or die "Can't output the output: $out-Merged.png because: $!\n";
	} else {
		open ($OFH, " | display - ") or die "Can't output to display.\n";
	}
	print $OFH $panel->png if $png==1;
	print $OFH $panel->svg if $svg==1;
	close $OFH;
}

sub help {
	die "Help: This script will take a 2 bams mapped @ different references and merge the reference & map at the merged reference. This is useful for merging LGT regions.
		--bam1=				bam1. Assumes position sorted & indexed.
		  --bam1_region=		<chr#:100-200>
		  --sort1=			<0|1> [0] 1= Position sort & index bam1.
		--bam2=				bam2. Assumes position sorted & indexed.
		  --bam2_region=		<chr#:100-200>
		  --sort2=			<0|1> [0] 1= Position sort & index bam2.
		--ref1=				Reference 1 fasta.
		  --ref1_region=		<chr#:100-200>
		--ref2=				Reference 2 fasta.
		  --ref2_region=		<chr#:100-200>
		--reads_list=			Path to a file with a list of desired reads to parse for. 1 read / line. 
		--M_only=			<0|1> [0] 1= When remapping to the merged reference, only keep M_* read pairs
		--MM_only=    			<0|1> [1] 1= When remapping to the merged reference, only keep M_M read pairs
		--merged_ref_name= 		Name for the new reference.  		[Merged]
		--n_num|n= 			Number of \"N's\" to insert inbetween merged references. [100]
		--png=				<0|1> [0] 1=Create a png img of the merged bam. 
		  --image_length=		Ajust the length of the png created.
		  --image_width=		Adjust the width of the png created. 
		  --pad_scale=			Pad white space around img. 
		  --fix_orientation=		<0|1> [1] 1= Try to determine how the references should be organized L-vs-R to make Mates face eachother.
							  0= Bam1 is on Left, Bam2 is on Right.
		--output_dir=			Directory for output. 				[bam1 dir]
		--output_prefix=		Prefix for the output fasta & bam. 	[bam1-merged-bam2]
		--help|h\n";
}

###################################

