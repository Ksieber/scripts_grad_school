#!/usr/bin/perl 
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use warnings;
no warnings 'uninitialized';
use strict;
use Data::Dumper;
use Carp;
use Bio::DB::Sam;
use Bio::Graphics;
use Bio::SeqIO;
use Bio::Util::DNA qw(:all);
use Bio::SeqFeature::Generic;
use GD;
use GD::SVG;
use POSIX;
use Math::NumberCruncher;
use Statistics::Distributions;
use Statistics::R;
use File::Basename;
use run_cmd;
use mk_dir;
use print_call;
use read_in_list;
use bwa;
use empty_chk;
use linecount;
use parse_flag;
use LGTSeek;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,             'input|bam1|i=s', 'bam2=s',         'ref1=s',            'ref2=s',            'bam1_region=s@',     'bam2_region=s@',  'draw_ref1_region=s@',
    'draw_ref2_region=s@', 'sort1=i',        'sort2=i',        'reads_list=s',      'fix_orientation=i', 'titrate_n_string=i', 'MM_only=i',       'anchor_bam1=i',
    'png=i',               'svg=i',          'draw_stdev|d=i', 'draw_both|B=i',     'picard_file|P=s',   'stdev|D=i',          'insert_size|I=i', 'image_length=i',
    'image_width=i',       'pad_scale=i',    'output_dir|o=s', 'output_prefix|p=s', 'merged_ref_name=s', 'dedup=i',            'jsd=i',           'threads|t=i',
    'Qsub|q=i',            'sub_mem=s',      'sub_name=s',     'number_of_reads=i', 'no_gal=i',          'help|?',             'min_cov=i',       'sub_mail=s',
    'max_window_size|w=i', 'more_info',
) or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} ) { &help; }    ## &help is at the end of the script
if ( $options{more_info} ) { &more_info; }       ## &more_info is at the end of the script

# Make sure we have the minimum inputs
if ( !$options{input} )       { die "Error: Must use: --bam1=<bam_aln_ref1> Please try agian.\n"; }
if ( !$options{picard_file} ) { die "Error: Must use: --picard_file=<LIB_insert_size_picard_file>. Please try again.\n"; }
if ( !$options{ref1} or !$options{ref2} ) { die "Error: Must use both: --ref1=<one_INT_ref.fa> --ref2=<other_INT_ref.fa>" }

# Ideally the user passes a bam with the reads already filtered for 1 INT
## The user can pass a bam with region(s) to pull reads from and try to create an INT model
## The user can specify a region to draw. This dictates both sides of the INT, but there may not be continious coverage. --draw_region overides --bam_region
## If there is no --bam_region and --draw_region, we will try to seperate the INT and calc each INT using &bam2regions_of_coverage.
## &bam2regions_of_coverage is quick naive look for potential regions with INT by identifying regions for bam1/ref1 coverage and bam2/ref2 with coverage.
## Based on the user passed or generated regions of intrest, we iterate over each region and analyze each pair of regions in-depth for INT and create a model foreach INT.

# Each list of regions = if (a file was passed and exists) ? (open the file to read in all the regions) : (else join all the different regions that can be split on "," )
my @draw_ref1_region_list;
if ( defined $options{draw_ref1_region} ) {
    @draw_ref1_region_list
        = ( -e "$options{draw_ref1_region}->[0]" )
        ? @{ &read_in_list( $options{draw_ref1_region}->[0] ) }
        : split( /,/, join( ',', @{ $options{draw_ref1_region} } ) );
    $options{draw_ref1_region} = join( ',', @draw_ref1_region_list );
}
my @draw_ref2_region_list;
if ( defined $options{draw_ref2_region} ) {
    @draw_ref2_region_list
        = ( -e "$options{draw_ref2_region}->[0]" )
        ? @{ &read_in_list( $options{draw_ref2_region}->[0] ) }
        : split( /,/, join( ',', @{ $options{draw_ref2_region} } ) );
    $options{draw_ref2_region} = join( ',', @draw_ref2_region_list );
}
my @bam1_region_list = undef;
if ( defined $options{bam1_region} ) {
    @bam1_region_list
        = ( -e "$options{bam1_region}->[0]" )
        ? @{ &read_in_list( $options{bam1_region}->[0] ) }
        : split( /,/, join( ',', @{ $options{bam1_region} } ) );
    $options{bam1_region} = join( ',', @bam1_region_list );
}
my @bam2_region_list = undef;
if ( defined $options{bam2_region} ) {
    @bam2_region_list
        = ( -e "$options{bam2_region}->[0]" )
        ? @{ &read_in_list( $options{bam2_region}->[0] ) }
        : split( /,/, join( ',', @{ $options{bam2_region} } ) );
    $options{bam2_region} = join( ',', @bam2_region_list );
}

# Submit the job to the SGE grid
if ( defined $options{Qsub} and $options{Qsub} == 1 ) {
    $options{sub_name} = defined $options{sub_name} ? $options{sub_name} : "mergeINTbams";
    Qsub_script( \%options );
}
print_call( \%options );

# Setup a few global default variables
my $JSD              = defined $options{jsd}              ? $options{jsd}              : 1;
my $MM_only          = defined $options{MM_only}          ? $options{MM_only}          : 1;
my $dedup            = defined $options{dedup}            ? $options{dedup}            : 1;
my $fix_orientation  = defined $options{fix_orientation}  ? $options{fix_orientation}  : 1;
my $titrate_n_string = defined $options{titrate_n_string} ? $options{titrate_n_string} : 1;
my $draw_nstring     = defined $options{draw_nstring}     ? $options{draw_nstring}     : 1;
my $draw_both        = defined $options{draw_both}        ? $options{draw_both}        : 1;
my $svg              = defined $options{draw_svg}         ? $options{draw_svg}         : 1;
my $png              = defined $options{draw_png}         ? $options{draw_png}         : "0";
my $min_coverage     = defined $options{min_cov}          ? $options{min_cov}          : 2;
my $threads          = defined $options{threads}          ? $options{threads}          : 3;

# Get filenames, paths, create output directory, and set output names.
my ( $fn1, $path1, $suff1 ) = fileparse( $options{input}, qr/\.[^\.]+/ );
my $out_dir = $options{output_dir} ? $options{output_dir} : $path1;
mk_dir("$out_dir");
my $output_prefix = $options{output_prefix} ? $options{output_prefix} : $fn1;
if ( $output_prefix =~ /(.+)\_psort$/ ) {
    $output_prefix = $1;    ## Remove a trailing _psort from input prefix
}
my $map_dir = "$out_dir/map_dir/";
mk_dir("$map_dir");
print_notebook( \%options );

## Map bam1 against ref2 if bam2 wasn't given.
print STDERR "======== Start: merge_2LGT_bams.pl ========\n";
if ( !$options{bam2} and defined $options{ref2} ) {
    print STDERR "======== BWA map bam1 at ref2 ========\n";
    my ( $fnR, $pathR, $suffR ) = fileparse( $options{ref2}, qr/\.[^\.]+/ );
    $options{bam2} = &bwa_aln(
        $options{input},
        $options{ref2},
        {   output_prefix => "$fn1\_$fnR",
            output_dir    => $map_dir,
            sort_index    => 1,
            cmd_log       => 1
        }
    );
}

# Coordinate sort the input bam if we need to
if ( defined $options{sort1} and $options{sort1} == 1 ) {
    print STDERR "======== Position sort bam1 ========\n";
    run_cmd("samtools sort -o $options{input} - | samtools view -b - > $map_dir/$fn1\_psort.bam");
    run_cmd("samtools index $map_dir/$fn1\_psort.bam");
    $options{input} = "$map_dir/$fn1\_psort.bam";
}
## Check to make sure we have coordinate sorted input.
else {
    my $header       = run_cmd("samtools view -H $options{input}");
    my @header_lines = split( /\n/, $header );
    my $hd_line      = $header_lines[0];
    if ( $hd_line !~ /SO\:coordinate/ ) {
        my ( $fn, $path, $suf ) = fileparse( $options{input}, ( ".nsrt.bam", ".srt.bam", qr/\.[^\.]+/ ) );
        run_cmd("samtools sort $options{input} $map_dir/$fn\.psrt");
        run_cmd("samtools index $map_dir/$fn\.psrt.bam");
        $options{input} = "$map_dir/$fn\.psrt.bam";
        $options{input} =~ s/\/{2,}/\//g;
    }
}

# Break apart bam1 into regions with coverage if not --bam_region was given
if ( !$options{bam1_region} ) {
    print STDERR "======== Parse regions of coverage for bam1 ========\n";
    @bam1_region_list = &bam2regions_of_coverage_v2(
        {   bam             => $options{input},
            ref             => $options{ref1},
            min_cov         => $min_coverage,
            max_window_size => $options{max_window_size}
        }
    );
}
print STDERR "======== split_bam1_regions:\n@bam1_region_list\n";

# Coordinate sort bam2
if ( defined $options{sort2} and $options{sort2} == 1 ) {
    print STDERR "======== Position sort bam2 ========\n";
    my ( $fn2, $path2, $suff2 ) = fileparse( $options{bam2}, qr/\.[^\.]+/ );
    run_cmd("samtools sort -o $options{bam2} - | samtools view -b - > $map_dir\/$fn2\_psort.bam");
    run_cmd("samtools index $map_dir\/$fn2\_psort.bam");
    $options{bam2} = "$map_dir\/$fn2\_psort.bam";
}

# Break apart bam2 into regions of coverage
if ( !$options{bam2_region} ) {
    print STDERR "======== Parse regions of coverage for bam2 ========\n";
    @bam2_region_list = &bam2regions_of_coverage_v2(
        {   bam             => $options{bam2},
            ref             => $options{ref2},
            min_cov         => $min_coverage,
            max_window_size => $options{max_window_size}
        }
    );
}
print STDERR "======== split_bam2_regions:\n@bam2_region_list\n";

# User can also pass a list of read-ids to use in calculating the INT
## Create a hash with the desired read-ids
my %reads;
if ( defined $options{reads_list} and -e $options{reads_list} ) {
    open( IN, "<", "$options{reads_list}" ) or confess "Can't open reads_list: $options{reads_list} because: $!\n";
    while (<IN>) {
        chomp;
        $reads{$_}++;
    }
    close IN;
}

### Process each of the INT regions
# The bam_regions we have are naive
## 1. Grab the reads from each bam_region
## 2. Merge bam1_region & bam2_region reads based on their support for the INT (PE reads, read1 maps ref1 only, read2 maps ref2 only)
## 3. Determine the model sequence and optimal distance between both sides of the INT
## 4. Create and map the merged.bam at the calculated model(s)
## 5. Draw the model(s)
foreach my $bam1_region (@bam1_region_list) {
    foreach my $bam2_region (@bam2_region_list) {
        my $dir;
        if ( $bam1_region and $bam2_region ) { $dir = "$out_dir/$bam1_region\_INT\_$bam2_region\/"; }
        elsif ($bam1_region) { $dir = "$out_dir/$bam1_region/"; }
        elsif ($bam2_region) { $dir = "$out_dir/$bam2_region/"; }
        else                 { $dir = $out_dir; }
        ## Cleanup the output_dir name
        $dir =~ s/:/_/g;
        $dir =~ s/\|/_/g;
        $dir =~ s/\/{2,}/\//g;
        mk_dir("$dir");

        print STDERR "========================================================================\n";
        print STDERR "======== Analyzing regions for INT: $bam1_region | $bam2_region ========\n";

        # 1. Grab important bam data from region1 of interest
        print STDERR "======== Pull bam1 data: $bam1_region ========\n";
        my $bam1_data = &pull_bam_data(
            {   bam          => $options{input},
                bam_region   => $bam1_region,
                draw_regions => \@draw_ref1_region_list,
            }
        );
        if ( $bam1_data->{'empty'} == 1 ) {
            print STDERR "======== Skipping empty bam1_region: $bam1_region\n";
            next;
        }

        my $bam2_data;
        if ( $options{bam2} ) {
            print STDERR "======== Pull bam2 data: $bam2_region ========\n";
            $bam2_data = &pull_bam_data(
                {   bam          => $options{bam2},
                    bam_region   => $bam2_region,
                    draw_regions => \@draw_ref2_region_list,
                }
            );
            if ( $bam2_data->{'empty'} == 1 ) {
                print STDERR "======== Skipping empty bam2_region: $bam2_region\n";
                next;
            }
        }

        # Create a bam with reads for this INT only
        print STDERR "======== Merge bams ========\n";
        my $merged_bam = &merge_bams(
            {   bam1_data  => $bam1_data,
                bam2_data  => $bam2_data,
                dedup      => $dedup,
                output_dir => "$dir/merged_bams/"
            }
        );
        if ( ( defined $merged_bam->{empty} and $merged_bam->{empty} == 1 ) or $merged_bam->{count} <= 2 ) {
            print STDERR "======== These regions do not have an INT: $bam1_region | $bam2_region ========\n";
            run_cmd("rm -rf $dir");
            next;
        }
        else {
            print STDERR "======== INT identified: bam1: $bam1_data->{bam_region} <==> bam2: $bam2_data->{bam_region} ========\n";
        }

        print STDERR "======== Create new INT-Ref ========\n";
        my $ref1;
        my $ref2;
        ## After this bam1 is always the "left"/5'/upstream bam
        if ( $fix_orientation == 1 ) {
            print STDERR "======== Fix INT-Ref Orientation ========\n";
            ( $bam1_data, $ref1, $bam2_data, $ref2 ) = &fix_orientation(
                {   merged_bam => $merged_bam,
                    bam1_data  => $bam1_data,
                    bam2_data  => $bam2_data,
                    ref1       => $options{ref1},
                    ref2       => $options{ref2},
                }
            );
        }

        # 2. Determine the sequences and size of the INT
        ## Calculate the optimal distance between the two sides of the INT using the AD & JSD
        print STDERR "======== Calculate INT-Ref sequence and distance ========\n";
        my ( $ref1_data_list, $ref2_data_list, @n_num_list ) = &optimize_refs(
            {   bam1_data        => $bam1_data,
                bam2_data        => $bam2_data,
                merged_bam       => $merged_bam,
                ref1             => $ref1,
                ref2             => $ref2,
                picard_file      => $options{picard_file},
                insert_size      => $options{insert_size},
                stdev            => $options{stdev},
                jsd              => $JSD,
                titrate_n_string => $titrate_n_string,
                MM_only          => $MM_only,
                output_dir       => $dir,
            }
        );

        # 4. Create and map the reads at the INT model
        print STDERR "======== Map INT reads at new-INT-Ref ========\n";
        my @bams_mapped_at_merged_refs_list = &merge_refs_and_map(
            {   bam1_data   => $bam1_data,
                bam2_data   => $bam2_data,
                ref1_region => $ref1_data_list,
                ref2_region => $ref2_data_list,
                merged_bam  => $merged_bam,
                min_cov     => $min_coverage,
                n_num       => \@n_num_list,
                output_dir  => $dir,
            }
        );

        # 5. Draw the INT model
        print STDERR "======== Draw new INT bam & Ref ========\n";
        foreach my $set (@bams_mapped_at_merged_refs_list) {
            ## Iterate to draw both STDEV & Normal img's if --draw_both
            if ( $draw_both == 1 ) {
                ## First draw the "normal" img
                &draw_img(
                    {   ref_data      => $set->{ref_data},
                        bam           => $set->{bam}->{file},
                        output_dir    => $set->{bam}->{dir},
                        output_prefix => $set->{bam}->{prefix},
                        png           => $png,
                        svg           => $svg,
                        draw_stdev    => "0",
                        draw_nstring  => $draw_nstring,
                    }
                );
                ## Second draw the color coded img
                &draw_img(
                    {   bam           => $set->{bam}->{file},
                        output_dir    => $set->{bam}->{dir},
                        output_prefix => "$set->{bam}->{prefix}\_stdev",
                        ref_data      => $set->{ref_data},
                        png           => $png,
                        svg           => $svg,
                        draw_stdev    => "1",
                        stdev         => $options{stdev},
                        insert_size   => $options{insert_size},
                        picard_file   => $options{picard_file},
                        draw_nstring  => $draw_nstring,
                    }
                );
            }
            elsif ( $options{png} or $options{svg} ) {
                &draw_img(
                    {   bam           => $set->{bam}->{file},
                        output_dir    => $set->{bam}->{dir},
                        output_prefix => $set->{bam}->{prefix},
                        ref_data      => $set->{ref_data},
                        png           => $png,
                        svg           => $svg,
                        draw_stdev    => $options{draw_stdev},
                        stdev         => $options{stdev},
                        insert_size   => $options{insert_size},
                        picard_file   => $options{picard_file},
                        draw_nstring  => $draw_nstring,
                    }
                );
            }
        }
        print STDERR "======== Finishd analysis for INT: $bam1_region | $bam2_region ========\n";
    }
}

run_cmd("rm -rf $map_dir");
print_complete( \%options );
## Done processing

###################################
########### Subroutines ###########
###################################

=head1

Title   : pull_bam_data
Usage   : my $bam_data = pull_bam_data($bam,$region);
Function: Extract the reads from a specific bam region and return them in an array. 
Args    : 
        bam             => /file/path/input.bam
        bam_region      => 'chr:100-200',
        draw_regions    => \@draw_ref#_region,

Returns : Returns a data structure with an array of the data as well as the data stored by hash.
        'file'          => /file/path/input.bam
        'id_hash'       => \%bam_data_hash,             ## $bam_data->{id_hash}->{$id}=$corresponding_bam_line_data
        'header'        => \@samtools_header
        'strand'        => Integer. Positive means more reads map to + strand.
        'bam_region'    => 'chr:100-200'
        'draw_regions'  => \@draw_ref#_region,          ## This links --draw_region to the object, allowing us to switch bam1 <-> bam2 later depending on the orientation and maintain the proper --draw_region <-> bam relationship
=cut

sub pull_bam_data {
    my $opts = shift;
    if ( &empty_chk( $opts->{bam} ) == 1 ) { confess "Error: The input: $opts->{bam} is empty.\n"; }

    ## Warn user if the specified region is empty
    if ( $opts->{bam_region} && run_cmd("samtools view -F0x4 $opts->{bam} \'$opts->{bam_region}\' | WC") == 0 ) {
        my $retval = { 'empty' => "1" };
        return $retval;
    }

    # Grab the samtools header, and remove the @PG line.
    my @header = `samtools view -H $opts->{bam}`
        or confess "Error: Couldn't get the samtools header from: $opts->{bam}\n";
    delete $header[-1];
    foreach my $header_lines (@header) {
        if ( $header_lines =~ /^\@HD/ ) { $header_lines =~ s/coordinate/queryname/; }
        if ( $header_lines =~ /^\@SQ\s+SN\:[A-Za-z0-9\|\.\_]+\:\d+\-\d+\s+LN\:\d+/ ) {
            confess "Error: This bam has a chromosome name that is invalid and will cause problems. Please remove or mask \":\\d+\-\\d+\" in the chr name: $header_lines\n";
        }
    }

    my $open_cmd
        = defined $opts->{bam_region}
        ? "samtools view -F0x4 $opts->{bam} \'$opts->{bam_region}\' |"
        : "samtools view -F0x4 $opts->{bam} |";

    my %bam_data;
    open( BAM, "$open_cmd" ) or confess "Error: Unable to open input bam: $opts->{bam} because: $!\n";
    while (<BAM>) {
        chomp;
        my @f = split;
        ## If a read list was passed, make sure this is a read we want
        my $read_id = $f[0];
        ## Add bam data to hash by id
        $bam_data{$read_id} = $_;
    }
    close BAM;

    my $retval = {
        'file'         => "$opts->{bam}",
        'id_hash'      => \%bam_data,
        'header'       => \@header,
        'empty'        => "0",
        'draw_regions' => $opts->{draw_regions},
        'bam_region'   => $opts->{bam_region},
    };
    return $retval;
}

=head1

Title   : merge_bams
Function: Create a bam with reads supporting an INT    
Usage   : my $merged_bam = &merge_bams(
            {   bam1_data  => $bam1_data,           ## pull_bam_data object
                bam2_data  => $bam2_data,           ## pull_bam_data object
                dedup      => $dedup,               ## <0|1> 1= dereplicate the data with picard
                output_dir => "$dir/merged_bams/"
            }
Returns : my $merged_bam = 
            {
                file        => $ret_bam,
                ids         => $merged_ids,        ## $hash{$ids}
                count       => $count,
                bam1_strand => $bam1_orientation,
                bam2_strand => $bam2_orientation,
            };

=cut

sub merge_bams {
    my $opts      = shift;
    my $bam1_data = $opts->{bam1_data};
    my $bam2_data = $opts->{bam2_data};
    mk_dir("$opts->{output_dir}");

    # Grab ID's present & mapping in both
    my $merged_ids = &merge_hash_ids( [ $bam1_data->{id_hash}, $bam2_data->{id_hash} ] );

    # Create a header for the tmp sam
    open( SAM, ">", "$opts->{output_dir}\/Merged-reads-only.sam" )
        or confess "Error: Couldn't open output merged.sam: $opts->{output_dir}/Merged-reads-only.sam";
    print SAM @{ $bam1_data->{header} };
    if ( $opts->{bam2_data} ) {
        foreach my $header_lines ( @{ $bam2_data->{header} } ) {
            if ( $header_lines !~ /^\@HD/ ) { print SAM $header_lines; }
        }
    }

    my $bam1_orientation = 0;
    my $bam2_orientation = 0;

    # Print the reads to the tmp mappig bam by using id's
    foreach my $id ( keys %{$merged_ids} ) {
        my $bam1_line = $bam1_data->{id_hash}->{$id};
        print SAM "$bam1_line\n";
        my $raw1_flag = ( split /\t/, $bam1_line )[1];
        my $readable1_flag = parse_flag($raw1_flag);
        if   ( $readable1_flag->{qrev} == 1 ) { $bam1_orientation--; }
        else                                  { $bam1_orientation++; }

        if ( $opts->{bam2_data} ) {
            my $bam2_line = $bam2_data->{id_hash}->{$id};
            print SAM "$bam2_line\n";
            my $raw2_flag = ( split /\t/, $bam2_line )[1];
            my $readable2_flag = parse_flag($raw2_flag);
            if   ( $readable2_flag->{qrev} == 1 ) { $bam2_orientation--; }
            else                                  { $bam2_orientation++; }
        }
    }

    close SAM;

    # Check to make sure we have reads in the merged bam.
    if ( &empty_chk("$opts->{output_dir}\/Merged-reads-only.sam") == 1 ) {
        return { empty => "1" };
    }

    # Convert the tmp mapping sam -> bam.
    run_cmd("samtools view -S $opts->{output_dir}\/Merged-reads-only.sam -bo $opts->{output_dir}\/Merged-reads-only.bam");
    run_cmd("rm $opts->{output_dir}\/Merged-reads-only.sam");
    my $ret_bam = "$opts->{output_dir}\/Merged-reads-only.bam";

    if ( $opts->{dedup} == 1 ) {
        my $LGTseq       = LGTSeek->new2();
        my $filtered_bam = $LGTseq->prinseqFilterBam(
            {   input_bam    => "$opts->{output_dir}\/Merged-reads-only.bam",
                output_dir   => "$opts->{output_dir}",
                rm_low_cmplx => "0",
                dedup        => "1",
            }
        );
        $ret_bam = $filtered_bam->{bam};

        run_cmd("samtools view $ret_bam | cut -f1 | sort -u > $opts->{output_dir}\/Post_dedup_good_ids.list");
        my $good_ids_hash = $LGTseq->_read_ids( { list => "$opts->{output_dir}\/Post_dedup_good_ids.list" } );
        run_cmd("rm $opts->{output_dir}\/Post_dedup_good_ids.list");

        foreach my $id ( keys %$merged_ids ) {
            if ( !$good_ids_hash->{$id} ) { delete $merged_ids->{$id}; }
        }
    }

    # Check to make sure we have reads in the merged bam.
    my $count = wc($ret_bam);

    return {
        file        => $ret_bam,
        ids         => $merged_ids,
        count       => $count,
        bam1_strand => $bam1_orientation,
        bam2_strand => $bam2_orientation,
    };
}

=head1

Title   : pull_ref_data
Usage   : my $ref_data_list = pull_ref_data($bam_data_in_array,@draw_ref_region_list);
Function: Extract the reference sequence from a specific fasta region for all ref_regions
Args    : 
        ref         =>  <file>              {Mandatory}
        ref_region  =>  <@array_ref>        {Optional} {Priority}
        bam_data    =>  <object>            {Mandatory}
        merged_ids  =>  <hash>              {Mandatory}
Returns : array_of_hashes{'seq' => ref_seq_string, 'range' => 'actual_chr:100-200'}

=cut

sub pull_ref_data {
    my $opts = shift;
    if ( !$opts->{ref} or !$opts->{bam_data} or !$opts->{merge_ids} ) {
        confess "Error: &pull_ref_data must have the following args passed: ref, bam_data, & merged_ids. Please try again.\n";
    }
    if ( &empty_chk( $opts->{ref} ) == 1 ) {
        confess "Error: The ref: $opts->{ref} is empty.\n";
    }

    my @ref_data;
    if ( defined $opts->{ref_region} ) {
        foreach my $region ( @{ $opts->{ref_region} } ) {
            my $ref_lower_range;
            my $ref_upper_range;
            my $ref_chr;

            if ( $opts->{ref_region} ) {
                $opts->{ref_region} =~ /^([\w\-\|\.]+)\:(\d+)\-(\d+)$/;
                ( $ref_chr, $ref_lower_range, $ref_upper_range ) = ( $1, $2, $3 );
            }

            # Pull the region from the fasta reference
            my $ref_db = Bio::DB::Fasta->new( $opts->{ref} );
            my $seq = ( $opts->{bam_data}->{rvcmplt} == 1 ) ? $ref_db->seq( $ref_chr, $ref_upper_range => $ref_lower_range ) : $ref_db->seq( $ref_chr, $ref_lower_range => $ref_upper_range );

            push(
                @ref_data,
                {   'seq'   => $seq,
                    'range' => "$ref_chr\:$ref_lower_range\-$ref_upper_range"
                }
            );

        }
    }
    else {
        my $ref_lower_range;
        my $ref_upper_range;
        my $ref_chr;
        my $bam_lower_range;
        my $bam_upper_range;
        my $bam_chr;

        my $first_read = 0;
        foreach my $read_id ( keys %{ $opts->{merged_ids} } ) {
            if ( $first_read == 0 ) {
                my $bam_line = $opts->{bam_data}->{id_hash}->{$read_id};
                my ( $chr, $position, $seq ) = ( split /\t/, $bam_line )[ 2, 3, 9 ];
                $bam_chr         = $chr;
                $bam_lower_range = $position;
                $bam_upper_range = $position + length($seq);
                $first_read++;
            }
            else {
                my $bam_line = $opts->{bam_data}->{id_hash}->{$read_id};
                my ( $chr, $position, $seq ) = ( split /\t/, $bam_line )[ 2, 3, 9 ];
                next if ( $chr ne $bam_chr );
                $bam_lower_range = $position if ( $position < $bam_lower_range );
                $bam_upper_range = ( $position + length($seq) )
                    if ( ( $position + length($seq) ) > $bam_upper_range );
            }
        }

        # Use reference region if it exist.
        my $actual_lower_range
            = defined $opts->{ref_region}
            ? $ref_lower_range
            : $bam_lower_range;
        my $actual_upper_range
            = defined $opts->{ref_region}
            ? $ref_upper_range
            : $bam_upper_range;
        my $actual_chr = defined $opts->{ref_region} ? $ref_chr : $bam_chr;

        # Finally, pull the region from the fasta reference
        my $ref_db = Bio::DB::Fasta->new( $opts->{ref} );
        my $seq
            = ( $opts->{bam_data}->{rvcmplt} == 1 ) ? $ref_db->seq( $actual_chr, $actual_upper_range => $actual_lower_range ) : $ref_db->seq( $actual_chr, $actual_lower_range => $actual_upper_range );

        push(
            @ref_data,
            {   'seq'   => $seq,
                'range' => "$actual_chr\:$actual_lower_range\-$actual_upper_range"
            }
        );
    }
    return \@ref_data;
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
    for ( my $i = 0; $i < $n_num; $i++ ) {
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
    foreach my $hash ( @{$ref_array_of_hashes} ) {
        foreach my $ids ( keys %{$hash} ) {
            $seen_ids{$ids}++;
        }
    }
    my %merged_ids;
    foreach my $id ( keys %seen_ids ) {
        if ( $seen_ids{$id} == $number_of_hashes ) {
            $merged_ids{$id}++;
        }
    }
    return \%merged_ids;
}

=head1

Title   : fix_orientation
Usage   : ( $bam1_data, $ref1, $bam2_data, $ref2 ) = &fix_orientation(
                        {   bam1_data => $bam1_data,
                            ref1      => $options{ref1},
                            bam2_data => $bam2_data,
                            ref2      => $options{ref2}};
Function:   Check the strand information in bam_Data and re-orders the 2 data sets if needed so bam1 = upstream.
Args    : 
        bam1_data   =>  $bam1_data
        bam2_data   =>  $bam2_data
        ref1        => ref1_file
        ref2        => ref2_file

Returns : Input in proper order according to strand


=cut

sub fix_orientation {

    my $opts = shift;
    my ( $bam1_data, $ref1, $bam2_data, $ref2 );

    ## Fix orientation based on majority of the sequencing reads
    if ( $opts->{merged_bam}->{bam1_strand} <= 0 && $opts->{merged_bam}->{bam2_strand} >= 0 ) {
        $bam1_data            = $opts->{bam2_data};
        $ref1                 = $opts->{ref2};
        $bam2_data            = $opts->{bam1_data};
        $ref2                 = $opts->{ref1};
        $bam1_data->{rvcmplt} = 0;
        $bam2_data->{rvcmplt} = 0;
    }
    elsif ( $opts->{merged_bam}->{bam1_strand} >= 0 && $opts->{merged_bam}->{bam2_strand} <= 0 ) {
        $bam1_data            = $opts->{bam1_data};
        $ref1                 = $opts->{ref1};
        $bam2_data            = $opts->{bam2_data};
        $ref2                 = $opts->{ref2};
        $bam1_data->{rvcmplt} = 0;
        $bam2_data->{rvcmplt} = 0;
    }
    elsif ( $opts->{merged_bam}->{bam1_strand} >= 0 && $opts->{merged_bam}->{bam2_strand} >= 0 ) {
        $bam1_data            = $opts->{bam2_data};
        $ref1                 = $opts->{ref2};
        $bam2_data            = $opts->{bam1_data};
        $ref2                 = $opts->{ref1};
        $bam1_data->{rvcmplt} = 0;
        $bam2_data->{rvcmplt} = 1;
    }
    elsif ( $opts->{merged_bam}->{bam1_strand} <= 0 && $opts->{merged_bam}->{bam2_strand} <= 0 ) {
        $bam1_data            = $opts->{bam2_data};
        $ref1                 = $opts->{ref2};
        $bam2_data            = $opts->{bam1_data};
        $ref2                 = $opts->{ref1};
        $bam1_data->{rvcmplt} = 1;
        $bam2_data->{rvcmplt} = 0;
    }
    return ( $bam1_data, $ref1, $bam2_data, $ref2 );
}

=head1

Title   : &optimize_refs
Function: Calculate the ideal distance between the two integration (INT) references (refs) based on insert size (i_size). 
Returns : A list of of reference positions and a # of bp between both references to draw to illustrate the optimized INT configuration.
Usage   : my ( $ref1_data_list, $ref2_data_list, @n_num_list ) = &optimize_refs(
            {   bam1_data         => $bam1_data,
                bam2_data         => $bam2_data,
                merged_bam        => $merged_bam,
                ref1              => $ref1,
                ref2              => $ref2,
                picard_file       => /path/to/LIB_picard_insert_size_metrics.txt
                output_dir        => /path/for/output/
            });

Args    : 
            bam1_data         => $bam1_data             ( &pull_bam_data object )
            bam2_data         => $bam2_data             ( &pull_bam_data object )
            merged_bam        => $merged_bam            ( &merge_Bams object )
            ref1              => /path/to/ref_1.fa
            ref2              => /path/to/ref_2.fa
            picard_file       => /path/to/picard_insert_size_metrics.txt for LIB
            output_dir        => /path/for/output,
            MM_only           => <0|1> 1= Only use reads pairs that have both pairs map to the merged refernce
            jsd               => <0|1> 1= Use Jensen-Shannon Distance calculations to determine distance between the 2 INT refs.
            titrate_n_string  => <0|1> 1= Titrated the LIB_stdev distances between the 2 refs for visualizing opti distance. 
            insert_size       => LIB insert size    (overrides picard file parsing)
            stdev             => LIB stdev          (overrides picard file parsing)

Workflow:
            1. Init LIB i_size counts
            2. Init INT i_size counts when INT refs are adjacent
                2A. Create a reference with zero bases between ( "N_0" ) the two sides of the integration refs
                2B. Map merged.bam @ N_0 ref we created (2A). Calc INT i_sizes.
                2C. Init INT i_size counts from the picard file (2B).
            3. Titrate the optimal distance between the two references using JSD and AD
            4. Return a list of the ref sequences and opti-N

Input object structure:
        pull_bam_data object:
            'file'              => /file/path/input.bam
            'id_hash'           => \%bam_data_hash,             ## $bam_data->{id_hash}->{$id}=$corresponding_bam_line_data
            'header'            => \@samtools_header
            'strand'            => Integer. Positive means more reads map to + strand.
            'bam_region'        => 'chr:100-200'
            'rvcmplt'           => <0|1>

        merge_Bams object:
            'file'              => /path/to/file.bam,
            'ids'               => $hash{$ids}
            'count'             => number of reads,
            'bam1_strand'       => #, > 0 = more reads map positive strand, < 0 more reads map to the reverse strand
            'bam2_strand'       => #

=cut

sub optimize_refs {
    my $opts = shift;
    if (   !$opts->{bam1_data}
        || !$opts->{bam2_data}
        || !$opts->{ref2}
        || !$opts->{ref1}
        || !$opts->{merged_bam} )
    {
        confess "Error: Must pass &optimize_refs the following opts: bam1_data, bam2_data, ref1, ref2, merged_ids\n";
    }

    # Local variables
    my $bam1_data           = $opts->{bam1_data};
    my $bam2_data           = $opts->{bam2_data};
    my $merged_bam          = $opts->{merged_bam};
    my $MM_only             = defined $opts->{MM_only} ? $opts->{MM_only} : "1";
    my $jsd                 = defined $opts->{jsd} ? "$opts->{jsd}" : "0";
    my $tmp_optimal_ref_dir = "$opts->{output_dir}/tmp_optimal_ref_dir/";
    mk_dir($tmp_optimal_ref_dir);

    # Data to return
    my @ret_ref1_data_list;
    my @ret_ref2_data_list;
    my @n_num_list;

    # 1. Initialize count of the LIB i_size population
    ## This is based on the picard file for the LIB mapped at the appropriate reference
    my %LIB_count;
    my $LIB_stdev;
    my $LIB_median_insert_size;
    open( PIC, "<", "$opts->{picard_file}" )
        or confess "Error: Unable to open picard_file for reading: $opts->{picard_file}\n";
    my $header_1 = 1;
    while (<PIC>) {
        chomp( my $picard_data_line = $_ );
        if ( $picard_data_line =~ /^MEDIAN_INSERT_SIZE/ ) {
            chomp( my $picard_insert_sizes = <PIC> );
            my ( $fr_median_i_size, $fr_abs_deviation, $fr_mean_INT_i_sizes, $fr_stdev ) = ( split /\t/, $picard_insert_sizes )[ 0, 1, 4, 5 ];
            $LIB_stdev              = defined $opts->{stdev}       ? $opts->{stdev}       : $fr_abs_deviation;
            $LIB_median_insert_size = defined $opts->{insert_size} ? $opts->{insert_size} : $fr_median_i_size;
        }
        elsif ( $picard_data_line =~ /^insert_size/ ) { $header_1 = 0; next; }
        elsif ( $header_1 == 1 )                      { next; }
        elsif ( $picard_data_line =~ /^\d+/ ) {
            my ( $i_size, $fr, $rf, $tandem ) = split( /\t/, $_ );
            $LIB_count{$i_size} = $fr;
        }
    }
    close PIC;

    # 2. Initilize a count of the INT insert-size population with 0 bp between the two references ( N_0 )
    ## 2A First, make a reference with adjacent consensus sequences for ref1 & ref2

    # 2A.1 Create the consensus for the reads for each bam and respective bam_region
    print STDERR "======== Calculating the consensus sequences for each side of the INT ========\n";
    ## Bam1
    my $half_threads = floor( $threads / 2 );
    open( my $OUT_1_vcf_fh,
        "|-", "samtools view - -u | samtools sort -@ $half_threads -O bam -T tmp_bam1_sort - | samtools mpileup -uAf $opts->{ref1} - | bcftools call -m -O z > $tmp_optimal_ref_dir/bam1_ref1.vcf.gz" )
        or die "Error: Unable to open a filehandle_1 to the VCF command.\n";
    print $OUT_1_vcf_fh @{ $bam1_data->{header} };

    ## Bam2
    open( my $OUT_2_vcf_fh,
        "|-", "samtools view - -u | samtools sort -@ $half_threads -O bam -T tmp_bam2_sort - | samtools mpileup -uAf $opts->{ref2} - | bcftools call -m -O z > $tmp_optimal_ref_dir/bam2_ref2.vcf.gz" )
        or die "Error: Unable to open a filehandle_2 to the VCF command.\n";
    print $OUT_2_vcf_fh @{ $bam2_data->{header} };

    # Local variables to capture position data of the reads from the merged data.
    ## Some reads may be removed between regions_of_coverage (not INT specific) to merged_bam (INT specific).
    ## Losing reads may have altered the exact region of the INT so we recalculate it here.
    ## Bam1
    my $bam1_chr;
    my $bam1_min;
    my $bam1_max;
    ## Bam2
    my $bam2_chr;
    my $bam2_min;
    my $bam2_max;

    # Print the bam data foreach read to the VCF while capturing position data
    foreach my $read_id ( keys %{ $merged_bam->{ids} } ) {
        my $bam1_line = $bam1_data->{id_hash}->{$read_id};
        my $bam2_line = $bam2_data->{id_hash}->{$read_id};

        print $OUT_1_vcf_fh "$bam1_line\n";
        print $OUT_2_vcf_fh "$bam2_line\n";

        my @bam1_split = split( /\t/, $bam1_line );
        my @bam2_split = split( /\t/, $bam2_line );

        if ( !$bam1_chr ) { $bam1_chr = $bam1_split[2]; }
        if ( !$bam2_chr ) { $bam2_chr = $bam2_split[2]; }

        my $bam1_read_id_5position = $bam1_split[3];
        my $bam2_read_id_5position = $bam2_split[3];

        my $bam1_read_id_3position = $bam1_split[3] + length( $bam1_split[9] ) - 1;    ## Subtract 1 b/c of zero based counting
        my $bam2_read_id_3position = $bam2_split[3] + length( $bam2_split[9] ) - 1;    ## Subtract 1 b/c of zero based counting

        if ( !$bam1_min ) { $bam1_min = $bam1_read_id_5position; }
        if ( !$bam2_min ) { $bam2_min = $bam2_read_id_5position; }

        if ( !$bam1_max ) { $bam1_max = $bam1_read_id_3position; }
        if ( !$bam2_max ) { $bam2_max = $bam2_read_id_3position; }

        if ( $bam1_read_id_5position < $bam1_min ) { $bam1_min = $bam1_read_id_5position; }
        if ( $bam2_read_id_5position < $bam2_min ) { $bam2_min = $bam2_read_id_5position; }

        if ( $bam1_read_id_3position >= $bam1_max ) {
            $bam1_max = $bam1_read_id_3position;
        }
        if ( $bam2_read_id_3position >= $bam2_max ) {
            $bam2_max = $bam2_read_id_3position;
        }
    }

    close $OUT_1_vcf_fh;
    close $OUT_2_vcf_fh;

    # Index the VCF file
    run_cmd("tabix $tmp_optimal_ref_dir/bam1_ref1.vcf.gz");
    run_cmd("tabix $tmp_optimal_ref_dir/bam2_ref2.vcf.gz");

    # Create the consensus sequences
    run_cmd("samtools faidx $opts->{ref1} \'$bam1_chr\:$bam1_min\-$bam1_max\' | bcftools consensus $tmp_optimal_ref_dir/bam1_ref1.vcf.gz > $tmp_optimal_ref_dir/bam1_ref1.fa");
    run_cmd("samtools faidx $opts->{ref2} \'$bam2_chr\:$bam2_min\-$bam2_max\' | bcftools consensus $tmp_optimal_ref_dir/bam2_ref2.vcf.gz > $tmp_optimal_ref_dir/bam2_ref2.fa");

    # Determine if we need to flip the orientation of the consensus sequence in order to have INT reads facing eachother. If we have to flip it, make an output-note of it.
    ## Bam1
    my $bam1_region_0 = ( $bam1_data->{rvcmplt} == 1 ) ? "$bam1_chr\:$bam1_max\-$bam1_min" : "$bam1_chr\:$bam1_min\-$bam1_max";
    if ( $bam1_data->{rvcmplt} == 1 ) {
        open( OUT, ">", "$opts->{output_dir}/REVERSE_COMPLEMENTED_$bam1_region_0.txt" )
            or confess "Error: Unable to open output: $opts->{output_dir}/REVERSE_COMPLEMENTED.txt\n";
        print OUT "REVERSE_COMPLEMENTED: $bam1_region_0\n";
        close OUT;
    }
    ## Bam2
    my $bam2_region_0 = ( $bam2_data->{rvcmplt} == 1 ) ? "$bam2_chr\:$bam2_max\-$bam2_min" : "$bam2_chr\:$bam2_min\-$bam2_max";
    if ( $bam2_data->{rvcmplt} == 1 ) {
        open( OUT, ">", "$opts->{output_dir}/REVERSE_COMPLEMENTED_$bam2_region_0.txt" )
            or confess "Error: Unable to open output: $opts->{output_dir}/REVERSE_COMPLEMENTED.txt\n";
        print OUT "REVERSE_COMPLEMENTED: $bam2_region_0\n";
        close OUT;
    }

    # Open Bio::SeqIO to parse the seq to create the adjacent reference
    ## Bam1
    my $consensus1_fh      = Bio::SeqIO->new( -format => 'Fasta', -file => "$tmp_optimal_ref_dir/bam1_ref1.fa" );
    my $ref1_consensus     = $consensus1_fh->next_seq();
    my $ref1_consensus_seq = ( $bam1_data->{rvcmplt} == 1 ) ? $ref1_consensus->revcom()->seq() : $ref1_consensus->seq();
    ## Bam2
    my $consensus2_fh      = Bio::SeqIO->new( -format => 'Fasta', -file => "$tmp_optimal_ref_dir/bam2_ref2.fa" );
    my $ref2_consensus     = $consensus2_fh->next_seq();
    my $ref2_consensus_seq = ( $bam2_data->{rvcmplt} == 1 ) ? $ref2_consensus->revcom()->seq() : $ref2_consensus->seq();

    # Now that we have the sequence and region for both INT references, add them to the list of data to return so that we can draw it later
    push(
        @ret_ref1_data_list,
        {   'seq'   => $ref1_consensus_seq,
            'range' => $bam1_region_0,
        }
    );
    push(
        @ret_ref2_data_list,
        {   'seq'   => $ref2_consensus_seq,
            'range' => $bam2_region_0,
        }
    );
    push( @n_num_list, "0" );

    # Create the adjacent reference for the INT with N=0. This will allow us to accurately calculate the i_size for the INT with N=0.
    print STDERR "======== Create a new reference with both side of the INT adjacent to eachother =========\n";
    my $adjacent_model_fa = "$tmp_optimal_ref_dir/adjacent_model_refs.fa";
    open( my $REF, ">", $adjacent_model_fa ) or die "Error: &optimize_refs unable to open output model reference: $adjacent_model_fa\n";
    ## Print fasta header
    print $REF ">adjacent_model_ref::$bam1_region_0\__$bam2_region_0\n";
    ## Print fasta sequence
    print $REF $ref1_consensus_seq . $ref2_consensus_seq . "\n";
    close $REF;
    $consensus1_fh->close();
    $consensus2_fh->close();

    # 2B. Map the merged bam at the merged_N0_reference
    print STDERR "======== BWA aln INT reads to the N_0 INT-Ref =========\n";
    ## bwa index merged_N0_reference
    run_cmd("bwa index $adjacent_model_fa");
    ## bwa align & use Picard to calculate the i_size for the merged INT reads
    my $adjacent_model_refs_bam = &bwa_aln(
        $opts->{merged_bam}->{file},
        $adjacent_model_fa,
        {   output_prefix  => "adjacent_model_refs",
            output_dir     => $tmp_optimal_ref_dir,
            insert_metrics => 1,
            MM_only        => $MM_only,
            cmd_log        => 1
        }
    );

    # 2C. Init INT_count by parsing the i_size data from the picard file for N_0
    print STDERR "======== Calculating the insert size for the INT reads aligned to the N_0 INT-Ref =========\n";
    my $adjacent_model_refs_insert_size_file = "$tmp_optimal_ref_dir/adjacent_model_refs\_std_insert.metrics";
    my @INT_i_sizes;
    open( my $int_N0_picard_fh, "<", "$adjacent_model_refs_insert_size_file" )
        or confess "Error: Unable to open picard_file for reading: $adjacent_model_refs_insert_size_file\n";
    ## Start reading the INT i_size data from the merged_N0_reference.bam picard file
    my $header_2 = 1;
    while (<$int_N0_picard_fh>) {
        chomp( my $picard_data_line = $_ );
        if ( $picard_data_line =~ /^insert_size/ ) { $header_2 = 0; next; }
        elsif ( $header_2 == 1 ) { next; }
        elsif ( $picard_data_line =~ /^\d+/ ) {
            my ( $i_size, $fr, $rf, $tandem ) = split( /\t/, $_ );
            for ( my $i = 1; $i <= $fr; $i++ ) {
                push( @INT_i_sizes, $i_size );
            }
        }
    }
    close $int_N0_picard_fh;

    # 3. Titrate the optimal distance between the two sides of the INT
    ## Open the output we will print the AD & JSD calculations to
    open( VAR, ">", "$opts->{output_dir}\/Variance_from_avg.txt" )
        or confess "Error: Unable to open output file: $opts->{output_dir}/Variance_from_avg.tx";
    ## Print Header
    printf VAR ( "%-20s%-20s", "N", "Variance_from_avg" );
    ## JSD calculation header
    if ( $jsd == 1 ) { printf VAR ( "%-20s%-20s%-20s", "JSD", 'ci_min', 'ci_max' ); }
    print VAR "\n";

    # Local hash variables to store the titration data
    my %AD_titration;     ##  AD{$N} = calc_AD_at_N
    my %jsd_titration;    ## JSD{$N} = calc_JSD_at_N

    # Titrate 0-100 bp between the INT refs, calculate the AD & JSD.
    print STDERR "======== Titrating the optimal distance between the consensus sequences =========\n";
    for ( my $N = 0; $N <= 100; $N++ ) {
        my @diff_N_list;    ## difference between the LIB_median_i_size and each INT_read_i_size w/ #_N's bp between the refs
        my %INT_count;      ## Same data structure as LIB_count
        foreach my $insert (@INT_i_sizes) {
            $INT_count{ ( $insert + $N ) }++;
            push( @diff_N_list, ( abs( $insert + $N - $LIB_median_insert_size ) ) );
        }
        ## Calculate the Average Difference
        my $avg_diff_N = Math::NumberCruncher::Mean( \@diff_N_list );
        $AD_titration{$N} = $avg_diff_N;
        printf VAR ( "%-20s%-20.3f", $N, $avg_diff_N );

        ## Calculate the Jensen-Shannon Distance
        if ( $jsd == 1 ) {
            ## Load functions into R to calculate the JSD
            my $R = Statistics::R->new( r_bin => '/usr/local/bin/R' );
            $R->run('require(boot)');
            $R->run(
                'calc_JSD <- function(inMatrix, pseudocount=0.0000001, ...) {
                                KLD <- function(x,y) { sum(x *log(x/y)) }
                                JSD <- function(x,y) { sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2)) }
                                matrixColSize <- length( colnames( inMatrix ) )
                                matrixRowSize <- length( rownames( inMatrix ) )
                                colnames <- colnames( inMatrix )
                                resultsMatrix <- matrix( 0, matrixColSize, matrixColSize )

                                inMatrix = apply( inMatrix, 1:2, function(x) ifelse ( x==0, pseudocount, x ) )
                                for ( i in 1:matrixColSize ) {
                                    for ( j in 1:matrixColSize ) {
                                        resultsMatrix[ i, j ] = JSD( as.vector( inMatrix[ , i ] ), as.vector( inMatrix[ , j ] ) )
                                    }
                                }
                                colnames -> colnames( resultsMatrix ) -> rownames( resultsMatrix )
                                as.dist( resultsMatrix ) -> resultsMatrix
                                attr( resultsMatrix, "method" ) <- "dist"
                                return( resultsMatrix )
                            }'
            );

            # $R->run(
            #     'calc_JSD_boot_fxn <- function (x_df, index) {
            #                 tmp_df <- data.frame(x_df[,1],x_df[index,2])
            #                 return( calc_JSD(tmp_df) )
            #                 }'
            # );
            # basic ci index
            $R->run(
                'calc_JSD_boot_fxn <- function (x_df, index) {
                            tmp_df <- data.frame(x_df[index,])
                            return( calc_JSD(tmp_df) )
                            }'
            );

            ## Initialize R data.frame with LIB & INT count of i_size in R
            $R->run('LIB_count = numeric()');
            $R->run('INT_count = numeric()');
            foreach my $key ( sort { $a <=> $b } keys %LIB_count ) {
                my $LIB_count_at_N_key = $LIB_count{$key};
                my $INT_count_at_N_key = defined $INT_count{$key} ? $INT_count{$key} : "0.0000001";    ## JSD can't have counts=0
                $R->run("LIB_count = c( LIB_count, $LIB_count_at_N_key )");
                $R->run("INT_count = c( INT_count, $INT_count_at_N_key )");
            }
            $R->run('counts=data.frame(LIB_count,INT_count)');
            ## Calculate the proportion of each i_size in LIB & INT
            $R->run('ct=prop.table(as.matrix(counts), margin=2)');

            ## Calculate the Jensen-Shannon Distance & parse output
            my $JSD_lines = $R->run('calc_JSD(ct)');
            my @JSD_split = split( /\n/, $JSD_lines );
            my $calc_JSD  = ( split /\s+/, $JSD_split[1] )[1];
            $jsd_titration{$N} = $calc_JSD;
            printf VAR ( "%-20.5f", $calc_JSD );

            # Calculate the JSD confidence interval for the model & parse output
            my $jsd_ci_lower;
            my $jsd_ci_upper;
            ## Bootstrap the INT population while keeping the LIB population_freq intact
            $R->run("JSD_boot <- boot(ct, calc_JSD_boot_fxn, R=1000, stype = \"i\", parallel=\"multicore\", ncpus=$threads)");
            my $JSdist_ci = $R->run('boot.ci(JSD_boot, type="norm")');
            my $ci_data_line = ( split /\n/, $JSdist_ci )[8];
            if ( defined $ci_data_line ) {
                $ci_data_line =~ /\s+\((.+)\,\s+(.+)\)/;
                $jsd_ci_lower = $1;
                $jsd_ci_upper = $2;
                printf VAR ( "%-20.4f%-20.4f", $jsd_ci_lower, $jsd_ci_upper );
            }
            else {
                $jsd_ci_lower = "NULL";
                $jsd_ci_upper = "NULL";
                printf VAR ( "%-20s%-20s", $jsd_ci_lower, $jsd_ci_upper );
            }

            ## Close R instance
            $R->stop();
        }
        print VAR "\n";
    }
    close VAR;

    # Determine the optimal distance between the two sides of the INT
    my $opti_AD_N = &find_key_with_min_hash_value( \%AD_titration );                              ## Optimal AD distance
    my $opti_JSD_N = ( $jsd == 1 ) ? &find_key_with_min_hash_value( \%jsd_titration ) : undef;    ## Optimal JSD distance

    # Print the optimal distance
    open( OPT, ">", "$opts->{output_dir}/Opti_dist.txt" )
        || confess "Error: Unable to open file to record optimal distance between references: $opts->{output_dir}/Opti_dist.txt\n";
    print OPT "Opti_Distance: $opti_AD_N";
    if ( $jsd == 1 ) { print OPT "\tJSD_Distance: $opti_JSD_N | JSD_value: $jsd_titration{$opti_JSD_N}"; }
    print OPT "\n";
    close OPT;

    # Graph the AD & JSD titration data
    my $R = Statistics::R->new( r_bin => '/usr/local/bin/R' );
    $R->run("Table=read.table(\"$opts->{output_dir}\/Variance_from_avg.txt\", header=T, row.names=1)");
    $R->run('df=data.frame(x=seq(0,length(Table[,2])-1), diff=Table[,1], jsd=Table[,2], lwr=Table[,3], upr=Table[,4])');
    if ( $jsd == 1 ) {
        ## Create the JSD plot
        $R->run("pdf(file=\"$opts->{output_dir}\/JSD_plot.pdf\")");
        $R->run('plot( jsd~x, data=df, ylim=range(c(df$lwr,df$upr)), cex=.1)');
        $R->run('with( df, polygon(c(x,rev(x)), c(lwr,rev(upr)), col="grey75", border=FALSE))');
        $R->run('matlines( df[,1], df[,c(-1,-2)], lwd=c(4,2,2), lty=1, col=c("black","red","red"))');
        $R->run("abline( h=$jsd_titration{$opti_JSD_N}, col=\"magenta\")");
        $R->run("abline( v=$opti_JSD_N, col=\"magenta\")");
        $R->run('dev.off()');
    }
    ## Create the AD plot
    $R->run("pdf(file=\"$opts->{output_dir}\/AD_plot.pdf\")");
    $R->run('plot(diff~x, data=df, ylim=range(df$diff), cex=.3, pch=19, cex.axis=0.6, cex.lab=0.6, font=2)');
    $R->run("abline( h=$AD_titration{$opti_AD_N}, col=\"magenta\")");
    $R->run("abline( v=$opti_AD_N, col=\"magenta\")");
    $R->run('dev.off()');

    # Close R instance
    $R->stop();

    if ( $opts->{titrate_n_string} == 1 ) {
        my @stdev_titration = ( 2, 1, .5, 0, -.5, -1 );
        foreach my $deviation (@stdev_titration) {
            my $step
                = ( $jsd == 1 )
                ? ( $opti_JSD_N + ( $LIB_stdev * $deviation ) )
                : ( $opti_AD_N + ( $LIB_stdev * $deviation ) );
            if ( $step >= 0 ) { push( @n_num_list, $step ); }
        }
    }

    print STDERR "======== Optimal reference calculated ========\n";
    run_cmd("rm -rf $tmp_optimal_ref_dir");
    return ( \@ret_ref1_data_list, \@ret_ref2_data_list, @n_num_list );
}

=head1

Title   : merge_refs_and_map
Usage   : my @bams_mapped_at_merged_refs_list = &merge_refs_and_map(
                                {   bam1_data   => $bam1_data,
                                    bam2_data   => $bam2_data,
                                    ref1_region => $ref1_data_list,
                                    ref2_region => $ref2_data_list,
                                    merged_bam  => $merged_bam,
                                    n_num       => \@n_num_list,
                                    output_dir  => $out_dir,
                                });
Function: Takes the sequences from both references, creates a new merged ref, and maps the input bams to the reference. 
Args    : 
        bam1_data   =>  $bam1_data
        bam2_data   =>  $bam2_data
        ref1_region =>  \@ref1_region_list
        ref2_region =>  \@ref2_region_list
        merged_ids  =>  $merged_ids
        min_cov     =>  $min_cov
        n_num       =>  \@n_num_list
        output_dir  =>  $output_dir (for ref_ranges.txt)
        M_only      => <0|1>
        MM_ony      => <0|1>

Returns : Array of hashes{'ref_data'}->{file || seq || cords},'bam'->{file || dir || prefix}}

=cut

sub merge_refs_and_map {
    my $opts = shift;

    my $M_only  = defined $opts->{M_only}  ? $opts->{M_only}  : "0";
    my $MM_only = defined $opts->{MM_only} ? $opts->{MM_only} : "0";
    my $merged_ids = $opts->{merged_bam}->{ids};

    my $ref1_data;
    my $ref2_data;
    my @ret_list;

    foreach my $ref1_region ( @{ $opts->{ref1_region} } ) {
        foreach my $ref2_region ( @{ $opts->{ref2_region} } ) {
            foreach my $n_num ( @{ $opts->{n_num} } ) {
                my $n_string   = create_n_string($n_num);
                my $out_suffix = "$ref1_region->{range}\_N$n_num\_$ref2_region->{range}";
                $out_suffix =~ s/:/_/g;
                $out_suffix =~ s/\|/_/g;
                my $out_dir = "$opts->{output_dir}\/$out_suffix";
                mk_dir("$out_dir");
                my $working_dir = "$out_dir\/merged_refs/";
                mk_dir("$working_dir");

                ## Build the new reference Sequence
                my $new_ref = $ref1_region->{seq} . $n_string . $ref2_region->{seq};

                ## Calculate the cordinates for the new reference
                my $cords = {
                    ref1_lower => 1,
                    ref1_upper => length( $ref1_region->{seq} ),
                    n_lower    => length( $ref1_region->{seq} ) + 1,
                    n_upper    => length( $ref1_region->{seq} ) + length($n_string),
                    ref2_lower => length( $ref1_region->{seq} ) + length($n_string) + 1,
                    ref2_upper => length( $ref1_region->{seq} ) + length($n_string) + length( $ref2_region->{seq} ),
                };

                ## Open text file to print the cordinates of the references and the cordinates of the new image.
                open( TXT, ">", "$working_dir/Ref_img_cords.txt" )
                    or confess "Error: can't open the output file for the ref & img cordinates: $working_dir/Ref_img_cords.txt because: $!";
                print TXT "Reference\tImage_cords\n";
                print TXT "$ref1_region->{range}\t$cords->{ref1_lower}\-$cords->{ref1_upper}\n";
                print TXT "N$n_num\:1-" . length($n_string) . "\t$cords->{n_lower}\-$cords->{n_upper}\n";
                print TXT "$ref2_region->{range}\t$cords->{ref2_lower}\-$cords->{ref2_upper}\n";
                close TXT;

                ## Print new reference Sequence.
                my $new_ref_name = defined $options{merged_ref_name} ? $options{merged_ref_name} : "Merged";
                open( FASTA, ">", "$working_dir\/Merged-refs\_$out_suffix.fa" )
                    or confess "Error: can't open the output file for the new merged fasta: $working_dir\/Merged-tmp\_$out_suffix.fa because: $!";
                print FASTA "\>$new_ref_name\n";    ## Print new "chr" header line for .fasta
                print FASTA $new_ref . "\n";
                close FASTA;
                run_cmd("bwa index $working_dir\/Merged-refs\_$out_suffix.fa");

                # Map the tmp bam @ the new reference.
                &bwa_aln(
                    $opts->{merged_bam}->{file},
                    "$working_dir\/Merged-refs\_$out_suffix.fa",
                    {   output_prefix => "Merged\-final\_$out_suffix",
                        output_dir    => $working_dir,
                        M_only        => $M_only,
                        MM_only       => $MM_only,
                        cmd_log       => 1
                    }
                );

                ## Rudimentry check for min coverage
                if ( ( run_cmd("samtools view $working_dir\/Merged\-final\_$out_suffix\.bam | WC") / 2 ) < $min_coverage ) {
                    run_cmd("rm -rf $out_dir");
                    next;
                }

                # Make a position sorted final output bam in the /working/dir/ to be used to draw the img
                run_cmd("samtools sort $working_dir\/Merged\-final\_$out_suffix\.bam $out_dir\/Merged\-final\_$out_suffix\-psrt");
                my $psrt_bam = "$out_dir\/Merged\-final\_$out_suffix\-psrt.bam";

                push(
                    @ret_list,
                    (   {   ref_data => {
                                file  => "$opts->{output_dir}\/Merged-refs\_$out_suffix.fa",
                                seq   => $new_ref,
                                cords => $cords,
                            },
                            bam => {
                                file   => $psrt_bam,
                                dir    => $out_dir,
                                prefix => "Merged-final\_$out_suffix",
                            }
                        }
                    )
                );
            }
        }
    }
    return @ret_list;
}

sub draw_img {
    my $opts       = shift;
    my $bam        = $opts->{bam};
    my $merged_ref = $opts->{ref_data};
    my $png        = defined $opts->{png} ? $opts->{png} : "0";
    my $svg        = defined $opts->{svg} ? $opts->{svg} : "0";
    my $draw_stdev = defined $opts->{draw_stdev} ? $opts->{draw_stdev} : "0";

    if ( &empty_chk($bam) == 1 ) {
        confess "Error: The Merged-final.bam: $bam is empty.\n";
    }
    #############################################################
    ## Setup OUTPUT
    my ( $fn, $path, $suff ) = fileparse( $bam, '.bam' );
    my $out_dir  = defined $opts->{output_dir}    ? $opts->{output_dir}    : $path;
    my $out_pref = defined $opts->{output_prefix} ? $opts->{output_prefix} : $fn;
    my $out      = "$out_dir\/$out_pref";
    #############################################################
    ## Setupt image boundries
    my $image_width  = defined $options{image_width}  ? $options{image_width}  : "1000";
    my $image_length = defined $options{image_length} ? $options{image_length} : "800";
    my $pad_scale    = defined $options{pad_scale}    ? $options{pad_scale}    : "0";
    #############################################################
    ## Calculate scale size ##
    my $scale_size = defined $options{scale} ? $options{scale} : length( $merged_ref->{seq} );
    if ( $image_width < $scale_size ) { $image_width = $scale_size + 20 }
    #############################################################
    # Calculate & setup drawing StDev.
    my $insert_size;
    my $stdev;
    if ( $draw_stdev == 1 ) {
        if ( defined $opts->{insert_size} && defined $opts->{stdev} ) {
            $insert_size = $opts->{insert_size};
            $stdev       = $opts->{stdev};
        }
        elsif ( defined $opts->{picard_file} ) {
            my @lines = `head -n 8 $opts->{picard_file}`;
            if ( $lines[6] !~ /^MEDIAN_INSERT_SIZE/ ) {
                confess "Error: The Picard file doesn't look right. Please fix it and try agian.\n";
            }
            ( $insert_size, $stdev ) = ( split /\t/, $lines[7] )[ 0, 1 ];
        }
        else {
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
    if ( $svg == 1 ) { push( @panel_options, ( -image_class => 'GD::SVG' ) ) }
    my $panel = Bio::Graphics::Panel->new(@panel_options);
    #############################################################
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
    ## Tracks for references
    ## Left Reference
    my $ref1_feature = Bio::SeqFeature::Generic->new(
        -start => $merged_ref->{cords}->{ref1_lower},
        -end   => $merged_ref->{cords}->{ref1_upper},
    );
    $panel->add_track(
        $ref1_feature,
        -glyph   => 'generic',
        -bgcolor => 'red',
    );

    # ## N-string "contig"
    if ( $opts->{draw_nstring} ) {
        my $nstring_feature = Bio::SeqFeature::Generic->new(
            -start => $merged_ref->{cords}->{n_lower},
            -end   => $merged_ref->{cords}->{n_upper},
        );
        $panel->add_track(
            $nstring_feature,
            -glyph   => 'crossbox',
            -fgcolor => 'black',
            -bgcolor => '128,128,128'
        );
    }
    ## Right Reference
    my $ref2_feature = Bio::SeqFeature::Generic->new(
        -start => $merged_ref->{cords}->{ref2_lower},
        -end   => $merged_ref->{cords}->{ref2_upper},
    );
    $panel->add_track(
        $ref2_feature,
        -glyph   => 'generic',
        -bgcolor => 'blue',
    );
    #############################################################
    ## Track for bam data
    ## Initialize the bam db
    my $sam = Bio::DB::Sam->new(
        -bam          => $bam,
        -fasta        => $merged_ref->{file},
        -expand_flags => "1",
        -autoindex    => "1",
    );
    ## Add feature by reads_list and region later
    my @pairs = sort( $sam->features( -type => 'read_pair' ) );
    for my $pair (@pairs) {
        my ( $first_mate, $second_mate ) = $pair->get_SeqFeatures;
        my $track;
        if ( $draw_stdev == 1 ) {
            my $read_insert_size = $pair->length;
            my $variance         = $insert_size - $read_insert_size;
            my $color;
            if ( abs($variance) == 0 ) {
                $color = 'white';
            }
            elsif ( abs($variance) <= ( .5 * $stdev ) ) {
                $color = '255,175,175' if $variance < 0;    ## Red
                $color = '200,255,200' if $variance > 0;    ## Green
            }
            elsif ( abs($variance) <= $stdev ) {
                $color = '255,20,20'   if $variance < 0;    ## Red
                $color = '115,255,115' if $variance > 0;    ## Green
            }
            elsif ( abs($variance) <= ( 2 * $stdev ) ) {
                $color = '185,0,0' if $variance <= 0;       ## Red
                $color = '0,165,0' if $variance > 0;        ## Green
            }
            elsif ( abs($variance) > ( 2 * $stdev ) ) {
                $color = '100,0,0' if $variance <= 0;       ## Red
                $color = '0,75,0'  if $variance > 0;        ## Green
            }
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
                -bgcolor   => 'purple',
            );
        }
        $track->add_feature( [ $first_mate, $second_mate ] );
    }

    #############################################################
    ## Print out the image
    my $OFH;
    my $output_suffix = $png == 1 ? ".png" : ".svg";
    if ( $opts->{output_dir} || $opts->{output_prefix} ) {
        open( $OFH, ">", "$out$output_suffix" )
            or die "Can't output the output: $out-Merged.png because: $!\n";
    }
    elsif ( $opts->{stdout} ) {
        $OFH = *STDOUT;
    }
    else {
        open( $OFH, " | display - " ) or die "Can't output to display.\n";
    }
    print $OFH $panel->png if $png == 1;
    print $OFH $panel->svg if $svg == 1;
    close $OFH;
}

sub find_key_with_min_hash_value {
    my $hash = shift;
    my ( $key,   @keys ) = keys %$hash;
    my ( $small, @vals ) = values %$hash;

    for ( 0 .. $#keys ) {
        if ( $vals[$_] < $small ) {
            $small = $vals[$_];
            $key   = $keys[$_];
        }
    }
    $key;
}

sub bam2regions_of_coverage_v2 {
    my $opts            = shift;
    my $max_window_size = defined( $opts->{max_window_size} ) ? $opts->{max_window_size} : length( run_cmd("samtools view -F0x4 $opts->{bam} | head -n 1 | cut -f10") ) * 3;
    my $window          = {};
    my @ret_region_list;

    my $header       = run_cmd("samtools view -H $opts->{bam}");
    my @header_lines = split( /\n/, $header );
    my $hd_line      = $header_lines[0];
    if ( $hd_line !~ /SO\:coordinate/ ) {
        my ( $fn, $path, $suf ) = fileparse( $opts->{bam}, ( ".nsrt.bam", ".srt.bam", qr/\.[^\.]+/ ) );
        run_cmd("samtools sort $opts->{bam} $map_dir/$fn\.psrt");
        run_cmd("samtools index $map_dir/$fn\.psrt.bam");
        $opts->{bam} = "$map_dir/$fn\.psrt.bam";
        $opts->{bam} =~ s/\/{2,}/\//g;
    }

    $window->{bam} = "$opts->{bam}";

    open( my $infh, "-|", "samtools mpileup -f $opts->{ref} -A $opts->{bam}" )
        or confess "ERROR: Can't open: samtools mpileup $opts->{ref} -A $opts->{bam}\n";
    while (<$infh>) {
        chomp;
        my ( $chr, $current_position, $cov ) = ( split /\t/, $_ )[ 0, 1, 3 ];
        if ( !$window->{chr} ) {
            ## Initialize new window
            $window->{chr}           = $chr;
            $window->{$chr}->{start} = $current_position;
            $window->{$chr}->{end}   = $current_position;
        }
        elsif ( $window->{chr} ne $chr ) {
            ## Add region to the list of regions to return
            @ret_region_list = &_bam2regions_of_coverage_add_region( $window, @ret_region_list, );

            # Delete old window & initialize new window
            $window->{chr}           = $chr;
            $window->{$chr}->{start} = $current_position;
            $window->{$chr}->{end}   = $current_position;
        }
        elsif ( $current_position < $window->{$chr}->{end} + $max_window_size ) {
            $window->{$chr}->{end} = $current_position;
        }
        elsif ( $current_position >= $window->{$chr}->{end} + $max_window_size ) {
            ## Add region to the list of regions to return
            @ret_region_list = &_bam2regions_of_coverage_add_region( $window, @ret_region_list );

            # Delete old window & initialize new window
            $window->{chr}           = $chr;
            $window->{$chr}->{start} = $current_position;
            $window->{$chr}->{end}   = $current_position;
        }
    }
    close $infh;

    # Check if we need to add region to the list of regions to return
    @ret_region_list = &_bam2regions_of_coverage_add_region( $window, @ret_region_list );

    return @ret_region_list;
}

sub _bam2regions_of_coverage_add_region {
    my ( $window, @ret_region_list ) = @_;

    my $chr          = $window->{chr};
    my $region       = "$chr\:$window->{$chr}->{start}\-$window->{$chr}->{end}";
    my $int_coverage = run_cmd("samtools view -F0x4 $window->{bam} \'$region\' | WC");
    if ( $int_coverage >= $min_coverage ) {
        push( @ret_region_list, $region );
    }
    return @ret_region_list;
}

sub help {
    die "Help: This script will take a bam and 2 different references to calculate the optimal distance between two consensus fragments supporting an integration between the two refs. 
    --------------------------------------------------------------------------------------------------------------------------------------------------------
    --bam1=                 bam1. Assumes coordinate sorted & indexed. If not, use --sort=1. (Mandatory)
      --bam1_region=        <chr#:100-200> Pull reads only from this region. (Recommended; without it all regions >= --min_cov will be inspected.)
      --sort1=              <0|1> [0] 1= Position sort & index bam1.
    --bam2=                 bam2. Assumes position sorted & indexed. If not, use --sort=1. (Optional)
      --bam2_region=        <chr#:100-200> Pull reads only from this region. (Recommended; without it all regions >= --min_cov will be inspected.)
      --sort2=              <0|1> [0] 1= Position sort & index bam2.
    --picard_file|P=        < /path/to/file.txt > Picard insert metrics file.
      --insert_size|I=      < # > overides data in the LIB picard file.
      --stdev|D=            < # > overides data in the LIB picard file.
    --min_cov=              < # > [2] Min reads coverage over an integration.
    --max_window_size|w     < # > [3*read-length] Max window *extension* for finding regions of coverage.
    --reads_list=           Path to a file with a list of desired reads to parse for. 1 read ID / line. 
    --------------------------------------------------------------------------------------------------------------------------------------------------------
    --ref1=                 Reference 1 fasta.      (Assumes bam1 is already mapped aginst ref1)
      --ref1_region=        <chr#:100-200> Use this reference range to map & draw reads against. ** More info below **
    --ref2=                 Reference 2 fasta.      (If no --bam2, and --ref2 is used, --bam1 will be mapped against --ref2)
      --ref2_region=        <chr#:100-200> Use this reference range to map & draw reads against.  ** More info below **
    --------------------------------------------------------------------------------------------------------------------------------------------------------
    --jsd=                <0|1> [1] 1= Calculate Jensen-Shannon Distance between the LIB and INT i_size distributions.
    --titrate_n_string=   <0|1> [1] 1= Draw optimal N distance between refs. +/- {2, 1.5, 1.0, 0.5} * \$LIB_STDEV
    --fix_orientation=    <0|1> [1] 1= Try to determine how the references should be organized L-vs-R to make Mates face eachother.
                                    0= Bam1 is on Left, Bam2 is on Right.
    --MM_only=            <0|1> [1] 1= When remapping to the merged reference, only keep M_M read pairs (Highly recommended=1)
    --merged_ref_name=    Name for the new reference.         [Merged]
    --------------------------------------------------------------------------------------------------------------------------------------------------------
    --svg=                  <0|1> [1] 1= Create a svg img of the merged bam.
    --png=                  <0|1> [0] 1= Create a png img of the merged bam.
      --image_length=       Ajust the length of the png created.
      --image_width=        Adjust the width of the png created. 
      --pad_scale=          Pad white space around img. 
    --------------------------------------------------------------------------------------------------------------------------------------------------------
    --draw_both|B=          <0|1> [1] 1= Draw both a \"normal\" & stdev color coded img.
    --draw_stdev|d=         <0|1> [0] 1= Color code the reads based on # of STDEV from the median insert size;
                                +/- STDEV * 0.5=Light Red/Green ; 1=Red/Green ; 2=Dark Red/Green
    --------------------------------------------------------------------------------------------------------------------------------------------------------                              
    --output_dir|o=         Directory for output.               [/options/bam1/dir/]
    --output_prefix|p=      Prefix for the output fasta & bam.  [bam1-merged-bam2]
    --threads|t=            < # > [3] # threads to use.
    --Qsub|q=               <0|1> [0] 1= Submit the job to grid. 
    --sub_mem=
    --sub_name=
    --sub_mail=
    --help|?                Show --options
    --more_info             Show description how the script works and more info on the --options (work in progress)
    --------------------------------------------------------------------------------------------------------------------------------------------------------\n";
}

sub more_info {
    die "\n\t--more_info (work in progress):
    \t\tThis script will take a bam and 2 references. Primarily this script will be used to calculate a model for the distance between the sides of an INT. 
    \t\t--bam1= Assumes this is position sorted. If no --bam2 is given, bam1 will be mapped at ref2.
    \t\tIdeally, bam1 has been filtered for high confidence reads supporting 1 INT.
    \t\tThe user can pass samtools style (chrZ:1-100) --bam_regions to pull reads from the region of the bam and try to create an INT model.
    \t\tThe user can specify a --ref_region to dictate both sides of the INT. These --ref_regions may not have coverage. --ref_region overides --bam_region.
    \t\tIf neither --bam_regions or --ref_region was used, the script will try break apart the --bam1 & --bam2 by coverage >1x and look for regions supporting an INT.
    \t\t--min_cov dictates how many reads must map to both sides of the INT. This does not imply that the INT has a spot with --min_cov mpileup coverage.
    \n"
}

###################################
