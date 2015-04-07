#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
use warnings;
no warnings 'uninitialized';
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
use POSIX;
use Math::NumberCruncher;
use Statistics::Distributions;
use Statistics::R;
use File::Basename;
use run_cmd;
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
    \%options,             'input|bam1=s',     'bam2=s',          'ref1=s',         'ref2=s',        'bam1_region=s@', 'bam2_region=s@',    'draw_ref1_region=s@',
    'draw_ref2_region=s@', 'sort1=i',          'sort2=i',         'reads_list=s',   'n_num|n=s@',    'draw_nstring=i', 'fix_orientation=i', 'titrate_refs=i',
    'titrate_n_string=i',  'M_only=i',         'MM_only=i',       'anchor_bam1=i',  'png=i',         'svg=i',          'draw_stdev|d=i',    'draw_both|B=i',
    'picard_file|P=s',     'stdev|D=i',        'insert_size|I=i', 'image_length=i', 'image_width=i', 'pad_scale=i',    'output_dir|o=s',    'output_prefix|p=s',
    'merged_ref_name=s',   'dedup=i',          'jsd=i',           'threads|t=i',    'Qsub|q=i',      'sub_mem=s',      'sub_name=s',        'number_of_reads=i',
    'split_bam1_cov=i',    'split_bam2_cov=i', 'no_gal=i',        'alternative=s',  'help|?',        'min_cov=i',      'sub_mail=s',
) or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} ) { &help; }    ## &help is @ the end of the script

if ( !$options{input} || !$options{ref2} ) {
    die "Error: You MUST pass ATLEAST the folowing arguements: --bam1=<lgt_bam_aln_ref1> --ref2=<some.fa> Please try agian.\n";
}
if ( !$options{ref1} ) { $options{ref1} = "/local/projects-t3/HLGT/references/hg19/hg19.fa"; }

if ( $options{jsd} && !$options{picard_file} ) {
    die "Error: You must pass --picard_file=<insert_size_picard_file> when you use --jsd=1. Please try again.\n";
}
if ( ( $options{draw_stdev} or $options{draw_nstring} ) and !$options{png} ) { $options{svg} = 1; }
if ( $options{titrate_refs} || $options{draw_stdev} || $options{titrate_n_string} || $options{jsd} ) {
    if ( !$options{stdev} || !$options{insert_size} ) {
        if ( -e $options{picard_file} ) {
            open( PIC, "<", "$options{picard_file}" )
                or die "Error: unable to open picard_file for reading: $options{picard_file}\n";
            my @header = <PIC> . <PIC> . <PIC> . <PIC> . <PIC> . <PIC> . <PIC>;
            chomp( my $picard_insert_sizes = <PIC> );
            my ( $fr_median_i_size, $fr_abs_deviation, $fr_mean_isize, $fr_stdev ) = ( split /\t/, $picard_insert_sizes )[ 0, 1, 4, 5 ];
            $options{stdev}       = $fr_abs_deviation;
            $options{insert_size} = $fr_median_i_size;
        }
        else {
            die "Error: You must pass --stdev=<#> or --picard_file=</path/to/file.txt> when using --titrate_refs, --titrate_n_string, --jsd, --draw_stdev. Please try again.\n";
        }
    }
}
$options{threads} = defined $options{threads} ? $options{threads} : "1";

## list of regions = if (a file was passed and exists) ? (open the file to read in all the regions)    : (else join all the different args breaking on ","if needed)
my @draw_ref1_region_list;
if ( $options{draw_ref1_region} ) {
    @draw_ref1_region_list
        = ( -e "$options{draw_ref1_region}->[0]" )
        ? @{ &read_in_list( $options{draw_ref1_region}->[0] ) }
        : split( /,/, join( ',', @{ $options{draw_ref1_region} } ) );
    $options{draw_ref1_region} = join( ',', @draw_ref1_region_list );
}

my @draw_ref2_region_list;
if ( $options{draw_ref2_region} ) {
    @draw_ref2_region_list
        = ( -e "$options{draw_ref2_region}->[0]" )
        ? @{ &read_in_list( $options{draw_ref2_region}->[0] ) }
        : split( /,/, join( ',', @{ $options{draw_ref2_region} } ) );
    $options{draw_ref2_region} = join( ',', @draw_ref2_region_list );
}

my @bam1_region_list = undef;
if ( $options{bam1_region} ) {
    @bam1_region_list
        = ( -e "$options{bam1_region}->[0]" )
        ? @{ &read_in_list( $options{bam1_region}->[0] ) }
        : split( /,/, join( ',', @{ $options{bam1_region} } ) );
    $options{bam1_region} = join( ',', @bam1_region_list );
}
my @bam2_region_list = undef;
if ( $options{bam2_region} ) {
    @bam2_region_list
        = ( -e "$options{bam2_region}->[0]" )
        ? @{ &read_in_list( $options{bam2_region}->[0] ) }
        : split( /,/, join( ',', @{ $options{bam2_region} } ) );
    $options{bam2_region} = join( ',', @bam2_region_list );
}

if ( $options{Qsub} ) {
    $options{sub_name} = defined $options{sub_name} ? $options{sub_name} : "mergeLGTbams";
    Qsub_script( \%options );
}
print_call( \%options );

## Setup a few Default values
# my $lgtseek = LGTSeek->new2(\%options);
my $fix_orientation = defined $options{fix_orientation} ? $options{fix_orientation} : "1";
my $draw_png        = defined $options{draw_png}        ? $options{draw_png}        : "0";
my $dedup           = defined $options{dedup}           ? $options{dedup}           : "1";

## Get filenames, paths, create output directory and set output names.
my ( $fn1, $path1, $suff1 ) = fileparse( $options{input}, qr/\.[^\.]+/ );
my $out_dir = $options{output_dir} ? $options{output_dir} : $path1;
run_cmd("mkdir -p $out_dir");
my $log = "$out_dir\/merge_2lgt_bams.log";
run_cmd( "touch $log", $log );
my $output_prefix = $options{output_prefix} ? $options{output_prefix} : $fn1;
if ( $output_prefix =~ /(.+)\_psort$/ ) {
    $output_prefix = $1;    ## Remove a trailing _psort from input prefix
}
my $map_dir = "$out_dir/map_dir/";
run_cmd( "mkdir -p $map_dir", $log );
my $map_log = "$map_dir/log.txt";
run_cmd( "touch $map_log", $log );
print_notebook( \%options );

my $optimize = ( $options{titrate_n_string} == 1 || $options{titrate_refs} == 1 || $options{jsd} == 1 ) ? "1" : "0";
if ( $optimize == 1 ) {
    $options{draw_nstring}    = 1;
    $options{fix_orientation} = 1;
    $options{MM_only}         = 1;
    $options{svg}             = 1;
    $options{draw_both}       = 1;
    $options{anchor_left}     = defined( $options{anchor_left} ) ? $options{anchor_left} : "1";
}
print STDERR "======== Start: merge_2LGT_bams.pl ========\n";
## Map bam1 against ref2 if bam2 wasn't given.
if ( !$options{bam2} && $options{ref2} ) {
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

my $sort1 = defined $options{sort1} ? $options{sort1} : "0";
my $sort2 = defined $options{sort2} ? $options{sort2} : "0";

# Cordinate sort we the option was passed
if ( $sort1 == 1 ) {
    print STDERR "======== Position sort bam1 ========\n";
    run_cmd( "samtools sort -o $options{input} - | samtools view -b - > $map_dir/$fn1\_psort.bam", $map_log );
    run_cmd( "samtools index $map_dir/$fn1\_psort.bam",                                            $map_log );
    $options{input} = "$map_dir/$fn1\_psort.bam";
}

# Checking to make sure we have cordinate sorted input.
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

if ( $options{split_bam1_cov} == 1 ) {
    print STDERR "======== Parse regions of coverage for bam1 ========\n";
    @bam1_region_list = &bam2regions_of_coverage_v2(
        {   bam     => $options{input},
            ref     => $options{ref1},
            min_cov => $options{min_cov},
        }
    );
    print STDERR "==== split_bam1_regions:\n@bam1_region_list\n";
}

if ( $sort2 == 1 ) {
    print STDERR "======== Position sort bam2 ========\n";
    my ( $fn2, $path2, $suff2 ) = fileparse( $options{bam2}, qr/\.[^\.]+/ );
    run_cmd( "samtools sort -o $options{bam2} - | samtools view -b - > $map_dir\/$fn2\_psort.bam", $map_log );
    run_cmd( "samtools index $map_dir\/$fn2\_psort.bam",                                           $map_log );
    $options{bam2} = "$map_dir\/$fn2\_psort.bam";
}
if ( $options{split_bam2_cov} == 1 ) {
    print STDERR "======== Parse regions of coverage for bam2 ========\n";
    @bam2_region_list = &bam2regions_of_coverage_v2(
        {   bam     => $options{bam2},
            ref     => $options{ref2},
            min_cov => $options{min_cov},
        }
    );
    print STDERR "==== split_bam2_regions:\n@bam2_region_list\n";
}

my %reads;    ## Hash of desired read id's.
if ( $options{reads_list} ) {
    open( IN, "<", "$options{reads_list}" ) or confess "Can't open reads_list: $options{reads_list} because: $!\n";
    while (<IN>) {
        chomp;
        $reads{$_}++;
    }
    close IN;
}

### Process files
foreach my $bam1_region (@bam1_region_list) {
    foreach my $bam2_region (@bam2_region_list) {
        my $dir;
        if ( $bam1_region && $bam2_region ) { $dir = "$out_dir/$bam1_region\_LGT\_$bam2_region\/"; }
        elsif ($bam1_region) { $dir = "$out_dir/$bam1_region/"; }
        elsif ($bam2_region) { $dir = "$out_dir/$bam2_region/"; }
        else                 { $dir = $out_dir; }
        $dir =~ s/:/_/g;
        $dir =~ s/\|/_/g;
        $dir =~ s/\/{2,}/\//g;

        run_cmd( "mkdir -p $dir", $log );

        # Grab important bam data from region of interest
        print STDERR "======== Pull bam1 data ========\n";
        my $bam1_data = &pull_bam_data(
            {   bam          => $options{input},
                bam_region   => $bam1_region,
                draw_regions => \@draw_ref1_region_list,
                log          => $map_log
            }
        );
        if ( $bam1_data->{'empty'} == 1 ) {
            print STDERR "*** Warning *** Skipping empty bam1_region: $bam1_region\n";
            next;
        }
        my $bam2_data;
        if ( $options{bam2} ) {
            print STDERR "======== Pull bam2 data ========\n";
            $bam2_data = &pull_bam_data(
                {   bam          => $options{bam2},
                    bam_region   => $bam2_region,
                    draw_regions => \@draw_ref2_region_list,
                    log          => $map_log
                }
            );
            if ( $bam2_data->{'empty'} == 1 ) {
                print STDERR "*** Warning *** Skipping empty bam2_region: $bam2_region\n";
                next;
            }
        }

        print STDERR "======== Merge bams ========\n";
        my $merged_bam = &merge_bams(
            {   bam1_data  => $bam1_data,
                bam2_data  => $bam2_data,
                dedup      => $dedup,
                output_dir => "$dir/merged_bams/"
            }
        );
        if ( ( defined $merged_bam->{empty} and $merged_bam->{empty} == 1 ) or $merged_bam->{count} <= 2 ) {
            print STDERR "*** Warning *** Merged_bam was empty, skipping region(s): $bam1_region | $bam2_region \n";
            run_cmd( "rm -rf $dir", $log );
            next;
        }
        print STDERR "======== Create new LGT-Ref ========\n";

        my $ref1;
        my $ref2;
        ## After this bam1 is always the "left" bam
        if ( $optimize == 1 || $options{fix_orientation} == 1 ) {
            print STDERR "======== Fix LGT-Ref Orientation ========\n";
            ( $bam1_data, $ref1, $bam2_data, $ref2 ) = &fix_orientation(
                {   merged_bam => $merged_bam,
                    bam1_data  => $bam1_data,
                    bam2_data  => $bam2_data,
                    ref1       => $options{ref1},
                    ref2       => $options{ref2},
                }
            );
        }

        ##
        my $ref1_data_list;
        my $ref2_data_list;
        ## Maybe set before hand?  = (defined $options{n_num} && @{$options{n_num}}>1 ) ? push(@n_num_list,split( /,/, join( ',', @{ $options{ref1_region} } ) ) : "0";
        my @n_num_list;
        if ( $optimize == 1 ) {
            print STDERR "======== Pull LGT-Ref Seq && Calculate Optimal LGT-Ref Spacing ========\n";
            ( $ref1_data_list, $ref2_data_list, @n_num_list ) = &optimize_refs(
                {   titrate_refs     => $options{titrate_refs},
                    titrate_n_string => $options{titrate_n_string},
                    bam1_data        => $bam1_data,
                    bam2_data        => $bam2_data,
                    merged_bam       => $merged_bam,
                    ref1             => $ref1,
                    ref2             => $ref2,
                    picard_file      => $options{picard_file},
                    insert_size      => $options{insert_size},
                    stdev            => $options{stdev},
                    jsd              => $options{jsd},
                    output_dir       => $dir,
                    anchor_bam1      => $options{anchor_bam1},
                    n_num            => $options{n_num},
                    log              => $log,
                }
            );
        }
        else {
            print STDERR "======== Pull LGT-Ref Seq ========\n";

            # Grab ref seq. for each
            $ref1_data_list = &pull_ref_data(
                {   ref        => $options{ref1},
                    ref_region => \@draw_ref1_region_list,
                    bam_data   => $bam1_data,
                    merged_ids => \$merged_bam->{ids},
                }
            );

            ## Make return multiple refs & take list of ref regions
            $ref2_data_list = &pull_ref_data(
                {   ref        => $options{ref2},
                    ref_region => \@draw_ref2_region_list,
                    bam_data   => $bam2_data,
                    merged_ids => \$merged_bam->{ids},
                }
            );
            @n_num_list
                = ( defined $options{n_num} && @{ $options{n_num} } > 1 )
                ? push( @n_num_list, split( /,/, join( ',', @{ $options{ref1_region} } ) ) )
                : "0";
        }

        ## Make return multiple refs & take list of ref regions
        ## Move optimize into here
        print STDERR "======== Map LGT-Merged-Bam at New-LGT-Ref ========\n";
        my @bams_mapped_at_merged_refs_list = &merge_refs_and_map(
            {   bam1_data   => $bam1_data,
                bam2_data   => $bam2_data,
                ref1_region => $ref1_data_list,
                ref2_region => $ref2_data_list,
                merged_bam  => $merged_bam,
                min_cov     => $options{min_cov},
                n_num       => \@n_num_list,
                output_dir  => $dir,
                log         => $log,
            }
        );
        print STDERR "======== Draw New Mapping LGT-bam at LGT-Ref ========\n";
        foreach my $set (@bams_mapped_at_merged_refs_list) {
            ## Iterate to draw both STDEV & Normal img's if --draw_both
            if ( $options{draw_both} ) {
                ## First draw the "normal" img
                &draw_img(
                    {   ref_data      => $set->{ref_data},
                        bam           => $set->{bam}->{file},
                        output_dir    => $set->{bam}->{dir},
                        output_prefix => $set->{bam}->{prefix},
                        png           => $options{png},
                        svg           => $options{svg},
                        draw_stdev    => "0",
                        draw_nstring  => $options{draw_nstring},
                    }
                );
                ## Second draw the color coded img
                &draw_img(
                    {   bam           => $set->{bam}->{file},
                        output_dir    => $set->{bam}->{dir},
                        output_prefix => "$set->{bam}->{prefix}\_stdev",
                        ref_data      => $set->{ref_data},
                        png           => $options{png},
                        svg           => $options{svg},
                        draw_stdev    => "1",
                        stdev         => $options{stdev},
                        insert_size   => $options{insert_size},
                        picard_file   => $options{picard_file},
                        draw_nstring  => $options{draw_nstring},
                    }
                );
            }
            elsif ( $options{png} || $options{svg} ) {
                &draw_img(
                    {   bam           => $set->{bam}->{file},
                        output_dir    => $set->{bam}->{dir},
                        output_prefix => $set->{bam}->{prefix},
                        ref_data      => $set->{ref_data},
                        png           => $options{png},
                        svg           => $options{svg},
                        draw_stdev    => $options{draw_stdev},
                        stdev         => $options{stdev},
                        insert_size   => $options{insert_size},
                        picard_file   => $options{picard_file},
                        draw_nstring  => $options{draw_nstring},
                    }
                );
            }
        }
    }
}

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
        log             => '/path/to/log.txt'
Returns : Returns a data structure with an array of the data as well as the data stored by hash.
        'hash'          => \%bam_data_hash,
        'header'        => \@samtools_header
        'strand'        => Integer. Positive means more reads map to + strand.
        'draw_regions'  => \@draw_ref#_region,
=cut

sub pull_bam_data {
    my $opts = shift;
    if ( &empty_chk( $opts->{bam} ) == 1 ) { confess "Error: The input: $opts->{bam} is empty.\n"; }

    ## Warn user if the specified region is empty
    if ( $opts->{bam_region} && run_cmd("samtools view -F0x4 $opts->{bam} \"$opts->{bam_region}\" | wc -l") == 0 ) {
        my $retval = { 'empty' => "1" };
        return $retval;
    }

    # Grab the samtools header, and remove the @PG line.
    my @header = `samtools view -H $opts->{bam} 2>>$opts->{log}`
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
        ? "samtools view -F0x4 $opts->{bam} \"$opts->{bam_region}\" 2>>$opts->{log} |"
        : "samtools view -F0x4 $opts->{bam} 2>>$opts->{log} |";

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
        'id_hash'      => \%bam_data,
        'header'       => \@header,
        'empty'        => "0",
        'draw_regions' => $opts->{draw_regions},
    };
    return $retval;
}

=head1

Title   : merge_bams
Usage   : 
Function:  
Args    : 
        
        
Returns : 

=cut

sub merge_bams {
    my $opts      = shift;
    my $bam1_data = $opts->{bam1_data};
    my $bam2_data = $opts->{bam2_data};
    run_cmd( "mkdir -p $opts->{output_dir}", $log );

    # print STDERR "bam1: " . Dumper($bam1_data) . "\n";
    # print STDERR "bam2: " . Dumper($bam2_data) . "\n";

    # Grab ID's present & mapping in both
    my $merged_ids = &merge_hash_ids( [ $bam1_data->{id_hash}, $bam2_data->{id_hash} ] );

    # Create a fake header for the tmp sam
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
    run_cmd( "samtools view -S $opts->{output_dir}\/Merged-reads-only.sam -bo $opts->{output_dir}\/Merged-reads-only.bam", $log );
    run_cmd( "rm $opts->{output_dir}\/Merged-reads-only.sam",                                                              $log );
    my $ret_bam = "$opts->{output_dir}\/Merged-reads-only.bam";

    if ( $opts->{dedup} == 1 ) {
        my $lgtseek      = LGTSeek->new2();
        my $filtered_bam = $lgtseek->prinseqFilterBam(
            {   input_bam    => "$opts->{output_dir}\/Merged-reads-only.bam",
                output_dir   => "$opts->{output_dir}",
                rm_low_cmplx => "0",
                dedup        => "1",
            }
        );
        $ret_bam = $filtered_bam->{bam};

        run_cmd("samtools view $ret_bam | cut -f1 | sort -u > $opts->{output_dir}\/Post_dedup_good_ids.list");
        my $good_ids_hash = $lgtseek->_read_ids( { list => "$opts->{output_dir}\/Post_dedup_good_ids.list" } );
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
        ref         =>  <file>              {Mandatory}{Priority}
        ref_region  =>  <@array_ref>        {Optional} {Priority}
        bam_data    =>  <object>            {Mandatory}
        merged_ids  =>  <hash>              {Mandatory}
Returns : array_of_hashes{'seq' => ref_seq_string, 'range' => 'actual_chr:100-200'}

=cut

sub pull_ref_data {
    my $opts = shift;
    if ( !$opts->{ref} || !$opts->{bam_data} || !$opts->{merge_ids} ) {
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
            my $seq = $ref_db->seq( $ref_chr, $ref_lower_range => $ref_upper_range );

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

        # Use reference region is it exist.
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
        my $seq = $ref_db->seq( $actual_chr, $actual_lower_range => $actual_upper_range );

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

    ## Fix orientation based on majority of the sequencing reads if the options is passed
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

Title   : optimize_refs
Usage   : ( $ref1_data_list, $ref2_data_list, @n_num_list ) = &optimize_refs(
                        {   titrate_refs     => $options{titrate_refs},
                            bam1_data         => $bam1_data,
                            bam2_data         => $bam2_data,
                            merged_bam        => $merged_bam,
                            ref1              => $ref1,
                            ref2              => $ref2,
                            insert_size       => $options{insert_size},
                            stdev             => $options{stdev},
                            output_dir        => $out_dir
                        )};

Function: Calculate the ideal distance between the two ref based on insertsize. Then titrate around ideal to illustrate.
Args    : 
            titrate_refs      => $options{titrate_refs},
            titrate_n_string  => $options{titrate_n_string},
            bam1_data         => $bam1_data,
            bam2_data         => $bam2_data,
            merged_bam        => $merged_bam,
            ref1              => $ref1,
            ref2              => $ref2,
            picard_file       => /path/to/picard_insert_size_metrics.txt
            insert_size       => $options{insert_size},
            stdev             => $options{stdev},
            number_of_reads   => $options{number_of_reads},
            output_dir        => $out_dir,
            anchor_bam1       => $options{anchor_bam1},
            n_num             => $options{n_num},
            log               => $log,

Returns : 

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

    my $bam1_data = $opts->{bam1_data};
    my $bam2_data = $opts->{bam2_data};
    my $ref1_data;
    my $ref2_data;

    my $bam1_chr;
    my $bam2_chr;
    my @r1_length_list;
    my @r2_length_list;
    my @bam1_position_list;
    my @bam2_position_list;
    my $jsd = defined $opts->{jsd} ? "$opts->{jsd}" : "0";

    foreach my $id ( keys %{ $opts->{merged_bam}->{ids} } ) {
        my $bam1_line = $bam1_data->{id_hash}->{$id};
        my $bam2_line = $bam2_data->{id_hash}->{$id};

        my @bam1_split = split( /\t/, $bam1_line );
        my @bam2_split = split( /\t/, $bam2_line );

        if ( !$bam1_chr ) { $bam1_chr = $bam1_split[2]; }
        if ( !$bam2_chr ) { $bam2_chr = $bam2_split[2]; }

        push( @r1_length_list, length( $bam1_split[9] ) );
        push( @r2_length_list, length( $bam2_split[9] ) );

        push( @bam1_position_list, $bam1_split[3] );
        push( @bam2_position_list, $bam2_split[3] );
    }

    my ( $bam1_max, $bam1_min ) = Math::NumberCruncher::Range( \@bam1_position_list );
    my ( $bam2_max, $bam2_min ) = Math::NumberCruncher::Range( \@bam2_position_list );

    my $r1_length = Math::NumberCruncher::Mean( \@r1_length_list );
    my $r2_length = Math::NumberCruncher::Mean( \@r2_length_list );

    my $bam1_max_adj = $bam1_max + ( $r1_length - 1 );
    my $bam2_max_adj = $bam2_max + ( $r2_length - 1 );

    my @Isize;
    foreach my $id ( keys %{ $opts->{merged_bam}->{ids} } ) {
        my $bam1_line     = $bam1_data->{id_hash}->{$id};
        my $bam2_line     = $bam2_data->{id_hash}->{$id};
        my $bam1_position = ( split /\t/, $bam1_line )[3];
        my $bam2_position = ( split /\t/, $bam2_line )[3];
        push( @Isize, ( 1 + $bam1_max_adj - $bam1_position ) + ( 1 + ( $bam2_position + ( $r2_length - 1 ) ) - $bam2_min ) );
    }

    my %pop_count;
    if ( $jsd == 1 ) {
        open( PIC, "<", "$opts->{picard_file}" )
            or confess "Error: Unable to open picard_file for reading: $options{picard_file}\n";
        my @header              = <PIC> . <PIC> . <PIC> . <PIC> . <PIC> . <PIC> . <PIC>;
        my @picard_insert_sizes = <PIC> . <PIC> . <PIC>;
        my $foobar              = <PIC> . <PIC> . <PIC>;
        while (<PIC>) {
            next if ( $_ !~ /^\d+/ );
            my ( $i_size, $fr, $rf, $tandem ) = split( /\t/, $_ );
            $pop_count{$i_size} = $fr;
        }
        close PIC;
    }

    open( VAR, ">", "$opts->{output_dir}\/Variance_from_avg.txt" )
        or confess "Error: Unable to open output file: $opts->{output_dir}/Variance_from_avg.tx";
    ## Print Header
    printf VAR ( "%-20s%-20s", "N", "Variance_from_avg" );
    ## JSD calculation header
    if ( $jsd == 1 ) { printf VAR ( "%-20s%-20s%-20s", "JSD", 'ci_min', 'ci_max' ); }
    print VAR "\n";

    my %tprob_titration;
    my %jsd_titration;
    my %variance_titration;

    for ( my $N = 0; $N <= 100; $N++ ) {
        my @adjusted_insert_size;    ## Read length + #N bp inbetween
        my @variance_N_list;         ## "variance" meaning the difference between each read's i_size w/ N's bp inbetween & the populations's i_size
        my %model_count;
        foreach my $insert (@Isize) {
            push( @adjusted_insert_size, ( $insert + $N ) );
            $model_count{ ( $insert + $N ) }++;
            push( @variance_N_list, ( abs( $insert + $N - $opts->{insert_size} ) ) );
        }
        ## Calculate average number of bases different than the population
        my $avg_N_variance = Math::NumberCruncher::Mean( \@variance_N_list );
        $variance_titration{$N} = $avg_N_variance;
        printf VAR ( "%-20s%-20.3f", "$N", "$avg_N_variance" );

        ## Calculate the Jensen-Shannon Divergence
        if ( $jsd == 1 ) {
            ## Load packages to calculate Jensen-Shannon Distance between populations
            my $R = Statistics::R->new( r_bin => '/usr/local/bin/R' );
            $R->run('require(boot)');
            $R->run(
                'calc_JSD <- function(inMatrix, pseudocount=0.0000001, ...) {
                                KLD <- function(x,y) { sum(x *log(x/y)) }
                                JSD<- function(x,y) { sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2)) }
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

            ## Initialize ref & model
            $R->run('ref_count   = numeric()');
            $R->run('model_count = numeric()');
            foreach my $key ( sort { $a <=> $b } keys %pop_count ) {
                my $pop_isize_count = $pop_count{$key};
                my $model_isize_count = defined $model_count{$key} ? $model_count{$key} : "0.0000001";
                ## Add the pop and model count of each insert size to the matrix by col.
                $R->run("ref_count   = c( ref_count,   $pop_isize_count   )");
                $R->run("model_count = c( model_count, $model_isize_count )");
            }
            $R->run('counts=data.frame(ref_count,model_count)');
            ## Calculate the proportion of each Isize in respective populations(ref & model)
            $R->run('ct=prop.table(as.matrix(counts), margin=2)');

            ## Calculate the Jensen-Shannon Distance & parse output
            my $JSD_lines = $R->run('calc_JSD(ct)');
            my @JSD_split = split( /\n/, $JSD_lines );
            my $calc_JSD  = ( split /\s+/, $JSD_split[1] )[1];
            $jsd_titration{$N} = $calc_JSD;
            printf VAR ( "%-20.5f", $calc_JSD );

            # Calculate the JSdist CI for the model & parse output
            my $jsd_ci_lower;
            my $jsd_ci_upper;
            ## Bootstrap the model population while keeping the reference population_freq intact
            $R->run(
                'calc_JSD_boot_fxn = function (x_df, index) {
                            tmp_df <- data.frame( x_df[,1], x_df[index,2] )
                            return ( calc_JSD(tmp_df) ) 
                            }'
            );
            $R->run("JSD_boot <- boot(ct, calc_JSD_boot_fxn, R=1000, parallel=\"multicore\", ncpus=$options{threads})");
            my $JSdist_ci = $R->run('boot.ci(JSD_boot, type="norm")');
            my $ci_data_line = ( split /\n/, $JSdist_ci )[8];
            if ( defined $ci_data_line ) {
                $ci_data_line =~ /\s+\((.+)\,\s+(.+)\)/;
                $jsd_ci_lower = $1;
                $jsd_ci_upper = $2;
                printf VAR ( "%-20.4f%-20.4f", $jsd_ci_lower, $jsd_ci_upper );    ##
            }
            else {
                $jsd_ci_lower = "NULL";
                $jsd_ci_upper = "NULL";
                printf VAR ( "%-20s%-20s", $jsd_ci_lower, $jsd_ci_upper );        ##
            }

            ## Close R instance
            $R->stop();
        }
        print VAR "\n";
    }
    close VAR;

    my $calculated_n = &find_key_with_min_hash_value( \%variance_titration );
    my $jsd_calc_n = ( $jsd == 1 ) ? &find_key_with_min_hash_value( \%jsd_titration ) : undef;

    if ( $jsd == 1 ) {
        my $R = Statistics::R->new( r_bin => '/usr/local/bin/R' );
        $R->run("Table=read.table(\"$opts->{output_dir}\/Variance_from_avg.txt\", header=T, row.names=1)");
        $R->run('df=data.frame(x=seq(0,length(Table[,2])-1), diff=Table[,1], jsd=Table[,2], lwr=Table[,3], upr=Table[,4])');
        ## Creat the JSD plot
        $R->run("pdf(file=\"$opts->{output_dir}\/JSD_plot.pdf\")");
        $R->run('plot( jsd~x, data=df, ylim=range(c(df$lwr,df$upr)), cex=.1)');
        $R->run('with( df, polygon(c(x,rev(x)), c(lwr,rev(upr)), col="grey75", border=FALSE))');
        $R->run('matlines( df[,1], df[,c(-1,-2)], lwd=c(4,2,2), lty=1, col=c("black","red","red"))');
        $R->run("abline( h=$jsd_titration{$jsd_calc_n}, col=\"magenta\")");
        $R->run("abline( v=$jsd_calc_n, col=\"magenta\")");
        $R->run('dev.off()');
        ## Creat the Diff plot
        $R->run("pdf(file=\"$opts->{output_dir}\/Diff_plot.pdf\")");
        $R->run('plot(diff~x, data=df, ylim=range(df$diff), cex=.3, pch=19, cex.axis=0.6, cex.lab=0.6, font=2)');
        $R->run("abline( h=$variance_titration{$calculated_n}, col=\"magenta\")");
        $R->run("abline( v=$calculated_n, col=\"magenta\")");
        $R->run('dev.off()');
        ## Close R instance
        $R->stop();
    }

    open( OPT, ">", "$opts->{output_dir}/Opti_dist.txt" )
        || confess "Error: Unable to open file to record optimal distance between references: $opts->{output_dir}/Opti_dist.txt\n";
    print OPT "Opti_Distance: $calculated_n";
    if ( $jsd == 1 ) { print OPT "\tJSD_Distance: $jsd_calc_n | JSD_value: $jsd_titration{$jsd_calc_n}"; }
    print OPT "\n";
    close OPT;

    my @n_num_list;
    my @calc_ref1_region_list;
    my @calc_ref2_region_list;

    my @stdev_titration = ( 2, 1, .5, 0, -.5, -1 );

    if ( $opts->{titrate_n_string} == 1 ) {
        foreach my $deviation (@stdev_titration) {
            my $step
                = ( $jsd == 1 )
                ? ( $jsd_calc_n + ( $opts->{stdev} * $deviation ) )
                : ( $calculated_n + ( $opts->{stdev} * $deviation ) );
            if ( $step >= 0 ) { push( @n_num_list, $step ); }
        }
        push( @calc_ref1_region_list, "$bam1_chr\:$bam1_min\-$bam1_max_adj" );
        push( @calc_ref2_region_list, "$bam2_chr\:$bam2_min\-$bam2_max_adj" );
    }
    elsif ( $opts->{titrate_refs} == 1 ) {
        if ( $opts->{anchor_bam1} == 1 ) {
            foreach my $deviation (@stdev_titration) {
                my $var = $opts->{insert_size} + ( $deviation * $opts->{stdev} );
                push( @calc_ref2_region_list, ( $bam2_chr . ":" . ( $bam2_min - $var ) . "-" . ($bam2_max_adj) ) );
            }
            push( @calc_ref1_region_list, "$bam1_chr\:$bam1_min\-$bam1_max_adj" );
        }
        else {
            foreach my $deviation (@stdev_titration) {
                my $var = $opts->{insert_size} + ( $deviation * $opts->{stdev} );
                push( @calc_ref1_region_list, ( $bam1_chr . ":" . ($bam1_min) . "-" . ( $bam1_max_adj + $var ) ) );
            }
            push( @calc_ref2_region_list, "$bam2_chr\:$bam2_min\-$bam2_max_adj" );
        }
        push( @n_num_list, "0" );
    }

    # Finally, pull the region from the fasta reference
    ## List to actual regions to return
    my @ref1_data_list;
    ## Allow --options{draw_ref#_region} to overide calculated regions.
    my @iter_ref1_region_list = ( defined @{ $bam1_data->{draw_regions} } ) ? @{ $bam1_data->{draw_regions} } : @calc_ref1_region_list;
    foreach my $region1 (@iter_ref1_region_list) {
        $region1 =~ /^([\w\-\|\.]+)\:(\d+)\-(\d+)$/;
        my ( $ref_chr, $ref_lower_range, $ref_upper_range ) = ( $1, $2, $3 );
        my $ref_db = Bio::DB::Fasta->new( $opts->{ref1} );
        my $seq
            = ( $bam1_data->{rvcmplt} == 1 )
            ? $ref_db->seq( $ref_chr, $ref_upper_range => $ref_lower_range )
            : $ref_db->seq( $ref_chr, $ref_lower_range => $ref_upper_range );
        my $region = ( $bam1_data->{rvcmplt} == 1 ) ? "$ref_chr\:$ref_upper_range\-$ref_lower_range" : "$ref_chr\:$ref_lower_range\-$ref_upper_range";
        if ( $bam1_data->{rvcmplt} == 1 ) {
            open( OUT, ">", "$opts->{output_dir}/REVERSE_COMPLEMENTED_$region1.txt" )
                or confess "Error: Unable to open output: $opts->{output_dir}/REVERSE_COMPLEMENTED.txt\n";
            print OUT "REVERSE_COMPLEMENTED: $region1\n";
            close OUT;
        }
        push(
            @ref1_data_list,
            {   'seq'   => $seq,
                'range' => $region,
            }
        );
    }

    ## List to actual regions to return
    my @ref2_data_list;
    ## Allow --options{draw_ref#_region} to overide calculated regions.
    my @iter_ref2_region_list = ( defined @{ $bam2_data->{draw_regions} } ) ? @{ $bam2_data->{draw_regions} } : @calc_ref2_region_list;
    foreach my $region2 (@iter_ref2_region_list) {
        $region2 =~ /^([\w\-\|\.]+)\:(\d+)\-(\d+)$/;
        my ( $ref_chr, $ref_lower_range, $ref_upper_range ) = ( $1, $2, $3 );
        my $ref_db = Bio::DB::Fasta->new( $opts->{ref2} );
        my $seq
            = ( $bam2_data->{rvcmplt} == 1 )
            ? $ref_db->seq( $ref_chr, $ref_upper_range => $ref_lower_range )
            : $ref_db->seq( $ref_chr, $ref_lower_range => $ref_upper_range );
        my $region = ( $bam2_data->{rvcmplt} == 1 ) ? "$ref_chr\:$ref_upper_range\-$ref_lower_range" : "$ref_chr\:$ref_lower_range\-$ref_upper_range";
        if ( $bam2_data->{rvcmplt} == 1 ) {
            open( OUT, ">", "$opts->{output_dir}/REVERSE_COMPLEMENTED_$region2.txt" )
                or confess "Error: Unable to open output: $opts->{output_dir}/REVERSE_COMPLEMENTED.txt\n";
            print OUT "REVERSE_COMPLEMENTED: $region2\n";
            close OUT;
        }
        push(
            @ref2_data_list,
            {   'seq'   => $seq,
                'range' => $region,
            }
        );
    }

    return ( \@ref1_data_list, \@ref2_data_list, @n_num_list );
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
                                    log         => $log,
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
        log         =>  "/path/to/log/file.txt"

Returns : Array of hashes{'ref_data'->{file || seq || cords},'bam'->{file || dir || prefix}}

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
                run_cmd( "mkdir -p $out_dir", $opts->{log} );

                my $working_dir = "$out_dir\/merged_refs/";
                run_cmd( "mkdir -p $working_dir", $opts->{log} );

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
                run_cmd( "bwa index $working_dir\/Merged-refs\_$out_suffix.fa 2>>$opts->{log}", $opts->{log} );

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
                if ( ( $opts->{min_cov} ) && ( run_cmd("samtools view $working_dir\/Merged\-final\_$out_suffix\.bam | wc -l ") / 2 ) < $opts->{min_cov} ) {
                    run_cmd("rm -rf $out_dir");
                    next;
                }

                # Make a position sorted final output bam in the /working/dir/ to be used to draw the img
                run_cmd( "samtools sort $working_dir\/Merged\-final\_$out_suffix\.bam $out_dir\/Merged\-final\_$out_suffix\-psrt", $log );
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
    my $max_window_size = defined( $opts->{max_window_size} ) ? $opts->{max_window_size} : "30";
    my $window          = {};
    my @ret_region_list;

    my $header       = run_cmd("samtools view -H $opts->{bam}");
    my @header_lines = split( /\n/, $header );
    my $hd_line      = $header_lines[0];
    if ( $hd_line !~ /SO\:coordinate/ ) {
        my ( $fn, $path, $suf ) = fileparse( $opts->{bam}, ( ".nsrt.bam", ".srt.bam", qr/\.[^\.]+/ ) );    ## FOOBAR
        run_cmd("samtools sort $opts->{bam} $map_dir/$fn\.psrt");
        run_cmd("samtools index $map_dir/$fn\.psrt.bam");
        $opts->{bam} = "$map_dir/$fn\.psrt.bam";
        $opts->{bam} =~ s/\/{2,}/\//g;
    }

    open( my $infh, "-|", "samtools mpileup -f $opts->{ref} -A $opts->{bam}" )
        or confess "ERROR: Can't open: samtools mpileup $opts->{ref} -A $opts->{bam}\n";
    while (<$infh>) {
        chomp;
        my ( $chr, $current_position, $cov ) = ( split /\t/, $_ )[ 0, 1, 3 ];
        if ( !$window->{chr} ) {

            # Initialize new window
            $window->{chr}           = $chr;
            $window->{$chr}->{start} = $current_position;
            $window->{$chr}->{end}   = $current_position;
            $window->{coverage} += $cov;
        }
        elsif ( $window->{chr} ne $chr ) {

            # Add region to the list of regions to return
            @ret_region_list = &_bam2regions_of_coverage_add_region( $window, @ret_region_list );

            # Delete old window
            foreach my $keys ( keys %$window ) { undef $window->{$keys}; }

            # Initialize new window
            $window->{chr}           = $chr;
            $window->{$chr}->{start} = $current_position;
            $window->{$chr}->{end}   = $current_position;
            $window->{coverage} += $cov;
        }
        elsif ( $current_position - $window->{$chr}->{end} <= $max_window_size ) {
            $window->{$chr}->{end} = $current_position;
            $window->{coverage} += $cov;
        }
        elsif ( $current_position - $window->{$chr}->{end} > $max_window_size ) {

            # Add region to the list of regions to return
            @ret_region_list = &_bam2regions_of_coverage_add_region( $window, @ret_region_list );

            # Delete old window
            foreach my $keys ( keys %$window ) { undef $window->{$keys}; }

            # Initialize new window
            $window->{chr}           = $chr;
            $window->{$chr}->{start} = $current_position;
            $window->{$chr}->{end}   = $current_position;
            $window->{coverage} += $cov;
        }
    }
    close $infh;

    # Check if we need to add region to the list of regions to return
    @ret_region_list = &_bam2regions_of_coverage_add_region( $window, @ret_region_list );

    return @ret_region_list;
}

sub _bam2regions_of_coverage_add_region {
    my ( $window, @ret_region_list ) = @_;
    my $min_cov = defined( $options{min_cov} ) ? $options{min_cov} : "1";

    my $chr          = $window->{chr};
    my $total_bp     = $window->{$chr}->{end} - $window->{$chr}->{start};
    my $avg_coverage = ceil( $window->{coverage} / $total_bp );
    if ( $avg_coverage >= $min_cov ) {
        push( @ret_region_list, "$chr\:$window->{$chr}->{start}\-$window->{$chr}->{end}" );
    }
    return @ret_region_list;
}

sub help {
    die
        "Help: This script will take a 2 bams mapped @ different references and merge the references & map at the merged reference. The final new bam can be drawn as a png or svg. This is useful for merging and illustrating LGT regions.
    --------------------------------------------------------------------------------------------------------------------------------------------------------
    --bam1=                 bam1. Assumes position sorted & indexed. If not, use --sort=1. (Mandatory)
      --bam1_region=        <chr#:100-200> Pull reads only from this region. (Highly recommended)
      --split_bam1_cov=     <0|1> [0] 1= Try to pull reads from regions of the bam with coverage.
      --sort1=              <0|1> [0] 1= Position sort & index bam1.
    --bam2=                 bam2. Assumes position sorted & indexed. If not, use --sort=1. (Optional)
      --bam2_region=        <chr#:100-200> Pull reads only from this region. (Highly recommended)
      --split_bam2_cov=     <0|1> [0] 1= Try to pull reads from regions of the bam with coverage.
      --sort2=              <0|1> [0] 1= Position sort & index bam2.
    --------------------------------------------------------------------------------------------------------------------------------------------------------
    --ref1=                 Reference 1 fasta.      (Assumes bam1 is already mapped aginst ref1)[ /local/projects-t3/HLGT/references/hg19/hg19.fa ]
      --ref1_region=        <chr#:100-200> Use this reference range to map & draw reads against. ** More info below **
    --ref2=                 Reference 2 fasta.      (If no --bam2, and --ref2 is used, --bam1 will be mapped against --ref2)
      --ref2_region=        <chr#:100-200> Use this reference range to map & draw reads against.  ** More info below **
    --------------------------------------------------------------------------------------------------------------------------------------------------------
    --titrate_refs=         <0|1> [0] 1= Calculate the ideal distance of sequence between reads. ** More Info on --titrate_refs below **
      --anchor_bam1=        <0|1> [0] 1= Keep upstream region the same while titrating downstream bam region. 0 = anchor Right. 
      --titrate_n_string=   <0|1> [0] 1= Calculate ideal distance of unknown sequence between reads.
      --jsd=                <0|1> [0] 1= Calculate Jensen-Shannon Divergence between the model and reference insert size distributions.
    --fix_orientation=      <0|1> [1] 1= Try to determine how the references should be organized L-vs-R to make Mates face eachother.
                                      0= Bam1 is on Left, Bam2 is on Right.
    --------------------------------------------------------------------------------------------------------------------------------------------------------
    --reads_list=           Path to a file with a list of desired reads to parse for. 1 read ID / line. 
    --M_only=               <0|1> [0] 1= When remapping to the merged reference, only keep M_* read pairs
    --MM_only=              <0|1> [1] 1= When remapping to the merged reference, only keep M_M read pairs (Highly recommended)
    --merged_ref_name=      Name for the new reference.         [Merged]
    --n_num|n=              Number of \"N's\" to insert inbetween merged references. [0] May also take comma delimited list or multiple entries. 
    --draw_nstring=         <0|1> [0] 1= Draw the n-string \"contig\".
    --------------------------------------------------------------------------------------------------------------------------------------------------------
    --png=                  <0|1> [0] 1= Create a png img of the merged bam.
    --svg=                  <0|1> [0] 1= Create a svg img of the merged bam.
      --image_length=       Ajust the length of the png created.
      --image_width=        Adjust the width of the png created. 
      --pad_scale=          Pad white space around img. 
    --------------------------------------------------------------------------------------------------------------------------------------------------------
    --draw_both|B=          <0|1> [0] 1= Draw both a \"normal\" & stdev color coded img.
    --draw_stdev|d=         <0|1> [0] 1= Color code the reads based on # of STDEV from the median insert size;
                                +/- STDEV * 0.5=Light Red/Green ; 1=Red/Green ; 2=Dark Red/Green
    --insert_size|I=        < # > 
    --stdev|D=              < # >
    --picard_file|P=        < /path/to/file.txt > Picard insert metrics file. Must be used if --stdev && --insert_size are not used or if jsd is calculated.
    --------------------------------------------------------------------------------------------------------------------------------------------------------                              
    --output_dir|o=         Directory for output.               [/options/bam1/dir/]
    --output_prefix|p=      Prefix for the output fasta & bam.  [bam1-merged-bam2]
    --stdout=               <0|1> [0] 1= Output goes to STDOUT. Either pipe it into a \"display\" ( | display - ) or redirect it to a new file ( > new.img)
    --help|?
    --------------------------------------------------------------------------------------------------------------------------------------------------------
     ** More info  ** 
    --ref#_region = May also accept a file with 1 region / line. May also take comma delimited regions or multiple --ref_region entries. If no region given, defaults to bam region. 
    Either --optimize option will calculate the ideal distance between references and draw +2, +1, +.5,0,-.5, -1 stdev's of read insert size around calculated ideal. 
    When using either --optimize option, you must pass --stdev & --insert_size.
     --------------------------------------------------------------------------------------------------------------------------------------------------------
    Example: perl merge_2lgt_bams.pl --bam1=bam_with_lgt_reads.bam --sort1=1 --ref1=hg19.fa --ref2=bacteria_with_lgt.fa --ref1_region=chr1:100-500 --bam1_region=chr1:200-400 --MM_only=1 --svg=1
    This will pull reads from the {bam_with_lgt_reads.bam} @ {chr1:200-400} ; map them against {bacteria_with_lgt.fa} ; create a new reference from {chr1:100-500} & 
    the region where the reads mapped in the {bacteria_with_lgt.fa}. The reads that are {mapped-mapped} to the new reference are then kept in the final bam and drawn as an {svg}. 
    --------------------------------------------------------------------------------------------------------------------------------------------------------\n";
}

###################################

__END__
            # $R->run('JSdist=vegdist( ct, method = "jensen-shannon", useShrinkage = TRUE )');
            # my $JSdist = $R->run('summary(JSdist)');
            # my $JSdist_split = ( split /\n/, $JSdist )[1];
            # my ( $JSmedian, $JSmean ) = ( split /\s+/, $JSdist_split )[ 2, 3 ];
            # printf VAR ( "\t%20s\t%20s", $JSmedian, $JSmean );
            # $jsd_titration{$N} = $JSmedian;

                        # Calculate the mean insert size CI for the model & parse output
            # $R->run('library(boot)');
            # $R->run('calc_mean_fxn=function(x,indices){ return( mean( x[indices] ) ) }');
            # $R->set( 'model_raw', \@Isize );
            # $R->run("model_boot_sampling = boot( model_raw, calc_mean_fxn, R=100, parallel=\"multicore\", ncpus=$options{threads})");
            # my $model_ci_data = $R->run('boot.ci(model_boot_sampling, type="norm")');
            # my @model_ci_data_split = split( /\n/, $model_ci_data );
            # my $ci_lower;
            # my $ci_upper;

            # if ( !$model_ci_data_split[8] ) {
            #     if ( $model_ci_data_split[1] eq 'NULL' ) {
            #         $ci_lower = "NULL";
            #         $ci_upper = "NULL";
            #     }
            # }
            # if ( $model_ci_data_split[8] =~ /\s+\((.+)\,\s+(.+)\)/ ) {
            #     $ci_lower = $1;
            #     $ci_upper = $2;
            # }
            # printf VAR ("\t%20s","$ci_lower -> $ci_upper");
