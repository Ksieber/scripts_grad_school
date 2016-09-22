#!/usr/bin/perl
my $VERSION = "1.01";
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use warnings;
use strict;
use Carp;
use run_cmd;
use mk_dir;
use print_call;
use mk_dir;
use empty_chk;
use File::Basename;
use POSIX;
use Statistics::R;
use Math::NumberCruncher;
use Data::Dumper;
use LGTSeek;
use linecount;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,   'pic=s',       'bam=s',    'output_dir|o=s', 'output_prefix|p=s', 'ref=s',     'subdirs=i', 'sort_mem=s', 'boot=i', 'reads=s',
    'iterate=i', 'threads|t=i', 'Qsub|q=i', 'help|?',         'sub_name=s',        'sub_mem=s', 'project=s', 'sub_mail=s', 'Qsub_iterate|Q=i',
) or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} ) { &help; }
if ( !$options{bam} ) { die "Error: Must pass input with --bam. Please try again.\n"; }
if ( defined $options{reads} and ( $options{reads} !~ /\d+/ and $options{reads} !~ /\d+\-\d+/ ) ) { die "Error: Your input doesn't look right. --reads=# or --reads=#-#.\n"; }

# Determine the range of # of reads to titrate over
my ( $reads_min, $reads_max ) = ( defined $options{reads} ) ? split( /\-/, join( "-", $options{reads} ) ) : ( 2, 100000 );
if ( !$reads_max and $reads_min =~ /\d+/ ) { $reads_max = $reads_min; }
elsif ( $reads_max !~ /\d+/ and $reads_min !~ /\d+/ ) { die "Unable to determine how many reads to downsample from the bam: $options{bam}.\n"; }
my @reads_number_to_process = ( $reads_min .. $reads_max );    # If I want to change it to be even # of reads only use: grep { $_ % 2 == 0 } (#..#)

# Determine output directory and prefix
my $lgtseq = LGTSeek->new2();
my ( $fn, $path, $suf )
    = ( defined $options{pic} and -e $options{pic} )
    ? fileparse( $options{pic}, $lgtseq->{'suffix_regex'} )
    : fileparse( $options{bam}, $lgtseq->{'bam_suffix_list'} );
my $output_dir = defined $options{output_dir} ? $options{output_dir} : "$path/tmp_jsd_titration";
mk_dir("$output_dir");
my $output_prefix = defined $options{output_prefix} ? $options{output_prefix} : $fn;

# Set a few other global varialbes
my $threads  = defined $options{threads} ? $options{threads} : "3";
my $Picard   = "/home/ksieber/bin/java -d64 -Xmx3G -Xms3G -XX:MaxDirectMemorySize=4500M -XX:ParallelGCThreads=$threads -jar $lgtseq->{Picard_jar}";    # -XX:+UseSerialGC
my $sort_mem = defined $options{sort_mem} ? $options{sort_mem} : "500M";
if ( defined $options{Qsub} or defined $options{Qsub_iterate} ) { $options{sub_mem} = defined $options{sub_mem} ? $options{sub_mem} : "8G"; }
my $iterate   = defined $options{iterate} ? $options{iterate} : "100";
my $Bootstrap = defined $options{boot}    ? $options{boot}    : "100";

## Qsub
if ( defined $options{Qsub} and $options{Qsub} == 1 ) {
    ## First, we don't want multiple jobs overwriting the input bam's picard file
    ## We will calculate it before submitting the jobs and assign each job the picard file to use
    my $ref_picard = ( defined $options{pic} and -e $options{pic} ) ? $options{pic} : &bam2picard( $options{bam}, $output_dir );    # Returns picard_file path
    $options{pic} = $ref_picard;
    wc( $options{bam} );                                                                                                            # count the input bam before other jobs start
    ## Set Qsub specific variables
    my $sub_name = defined $options{sub_name} ? $options{sub_name} : "JSDcalc";

    # Launch 1 Job / desired_#_of_reads
    foreach my $number_of_reads (@reads_number_to_process) {
        my $sub_mail
            = ( ( defined $lgtseq->{sub_mail} ) and ( $lgtseq->{sub_mail} == 1 and $number_of_reads == $reads_number_to_process[-1] ) )
            ? 1
            : 0;    # Only set the last job to send a mail. It may not finish last, but it should be close.
        my $cmd = "$^X $0";
        $cmd = $cmd . " --reads=$number_of_reads";
        foreach my $key ( keys %options ) {
            next if ( $key eq 'Qsub' );
            next if ( $key eq 'reads' );      # We are manually setting the number of reads above (line 67), so we don't want to have --options{reads} being passed
            next if ( $key eq 'sub_mail' );
            if ( defined $options{$key} ) { $cmd = $cmd . " --$key=$options{$key}" }
        }
        Qsub(
            {   cmd      => $cmd,
                wd       => $output_dir,
                threads  => $threads,
                sub_name => $sub_name,
                sub_mem  => $options{sub_mem},
                sub_mail => $sub_mail,
                project  => $lgtseq->{project},
            }
        );
        sleep 2;
    }
    die "+++ Finished submitting the jobs to the grid. +++\n";
}
elsif ( defined $options{Qsub_iterate} and $options{Qsub_iterate} == 1 ) {
    ## Set a few Qsub specific variables
    $options{sub_name} = defined $options{sub_name} ? $options{sub_name} : "JSDcalc";
    $options{sub_mem} = $lgtseq->{sub_mem};
    Qsub_script( \%options );
}

## Start the actual analysis
print_call( \%options, "## $0\-Version:$VERSION" );
print_notebook( \%options, "## $0\-Version:$VERSION" );
## Get input bam insert size data
my $ref_picard = ( defined $options{pic} and -e $options{pic} ) ? $options{pic} : &bam2picard( $options{bam} );    # Returns picard_file path
my $ref_pop = &pop_insert_size_from_picard_file($ref_picard);                                                      # Returns Hash_ref->{insert_size}=count

# Open output for JSD calcs.
my $print_header;
if ( !-e "$output_dir/$output_prefix\_jsd_titration.txt" ) { $print_header = 1; }                                  # Determine if we need to print a header or if it has already been done
open( my $OUT, ">", "$output_dir/$output_prefix\_jsd_titration\_$options{reads}.txt" ) or die "Error: Unable to open output to: $output_dir/$output_prefix\_jsd_titration.txt\n";
$| = 1;
select $OUT;
if ( defined $print_header and $print_header == 1 ) { print $OUT "#_Reads\tJSD\tCI_min\tCI_max\n"; }

# Foreach # of reads in the range of desired [--reads]
foreach my $desired_number_of_reads (@reads_number_to_process) {
    my $downsampled_bams_dir = "$output_dir/downsampled_bams\_$desired_number_of_reads/";
    ## We downsample each # of [--reads] for 1.. <= [--iterate] times
    for ( my $n = 1; $n <= $iterate; $n++ ) {
        mk_dir($downsampled_bams_dir);
        my $model_bam_with_number_of_reads = &downsample_bam( $desired_number_of_reads, $downsampled_bams_dir, $options{bam} );
        my $actual_number_of_reads         = wc($model_bam_with_number_of_reads);
        my $model_picard                   = &bam2picard($model_bam_with_number_of_reads);
        my $model_pop                      = &pop_insert_size_from_picard_file($model_picard);
        my ( $jsd, $jsd_ci_min, $jsd_ci_max ) = &calc_jsd( $model_pop, $ref_pop );    ## Ref_pop must be second for bootstrapping CI to work properly, #version 1.01
        print $OUT "$actual_number_of_reads\t$jsd\t$jsd_ci_min\t$jsd_ci_max\n";
        run_cmd("rm -rf $downsampled_bams_dir");
    }
}

print_complete( \%options );
## COMPLETE

#################
## Subroutines ##
#################

sub bam2picard {
    ## Grab ARGS
    my $bam = shift;
    my $dir = shift;
    ## Make sure we got a bam at least
    if ( !$bam )    { die "Error: &bam2picard didn't get a bam file.\n"; }
    if ( !-e $bam ) { die "Error: &bam2picard input bam doesn't exist: $bam\n"; }
    ## Set output_dir
    my ( $out_fn, $in_path, $out_suf ) = fileparse( $bam, $lgtseq->{bam_suffix_list} );
    my $out_dir = defined $dir ? $dir : $in_path;
    ## Check if we need to position sort the bam
    my $header = run_cmd("samtools view -H $bam | head -n 1");
    if ( $header !~ /coordinate/ ) {
        run_cmd("samtools sort -@ $threads -m $sort_mem $bam $out_dir/$out_fn\.foobar");
        $bam = "$out_dir/$out_fn\.foobar.bam";
    }
    ## Calculate the picard insert metrics
CALC_INSERT_SIZES: run_cmd("$Picard CollectInsertSizeMetrics AS=true I=$bam H=/dev/null O=$out_dir\/$out_fn\_insert.metrics M=0 VALIDATION_STRINGENCY=SILENT 2>>$out_dir/picard_stderr.log");
    if ( !-e "$out_dir\/$out_fn\_insert.metrics" or wc("$out_dir\/$out_fn\_insert.metrics") == 0 ) { goto CALC_INSERT_SIZES; }
    ## Return the path to the picard insert metrics file
    return "$out_dir/$out_fn\_insert.metrics";
}

sub pop_insert_size_from_picard_file {
    my $picard_file = shift;

    my %reference;
    open( PIC, "<", "$picard_file" ) or confess "Error: Unable to open the picard_file: $picard_file\n";
    my @header              = <PIC> . <PIC> . <PIC> . <PIC> . <PIC> . <PIC> . <PIC>;
    my @picard_insert_sizes = <PIC> . <PIC> . <PIC>;
    my $foobar              = <PIC> . <PIC> . <PIC>;
    while (<PIC>) {
        next if ( $_ !~ /^\d+/ );
        my ( $i_size, $fr, $rf, $tandem ) = split( /\t/, $_ );
        $reference{$i_size} = $fr;
    }
    close PIC;
    return \%reference;
}

sub calc_jsd {
    my $pop1 = shift;
    my $pop2 = shift;
    ## Intialize R and JSD function
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
    ## Initialize data sets in R
    $R->run('pop_1_count   = numeric()');
    $R->run('pop_2_count = numeric()');
    if ( scalar( keys %$pop1 ) > scalar( keys %$pop2 ) ) {
        foreach my $insert ( sort { $a <=> $b } keys %$pop1 ) {
            my $pop_1_isize_count = $pop1->{$insert};
            my $pop_2_isize_count = defined $pop2->{$insert} ? $pop2->{$insert} : "0.000000000001";
            ## Add the pop and model count of each insert size to the matrix by col.
            $R->run("pop_1_count = c( pop_1_count, $pop_1_isize_count   )");
            $R->run("pop_2_count = c( pop_2_count, $pop_2_isize_count )");
        }
    }
    else {
        foreach my $insert ( sort { $a <=> $b } keys %$pop2 ) {
            my $pop_1_isize_count = defined $pop1->{$insert} ? $pop1->{$insert} : "0.000000000001";
            my $pop_2_isize_count = defined $pop2->{$insert} ? $pop2->{$insert} : "0.000000000001";
            ## Add the pop and model count of each insert size to the matrix by col.
            $R->run("pop_1_count = c( pop_1_count, $pop_1_isize_count   )");
            $R->run("pop_2_count = c( pop_2_count, $pop_2_isize_count )");
        }
    }

    $R->run('counts=data.frame(pop_1_count,pop_2_count)');
    ## Calculate the proportion of each Isize in each populations
    $R->run('ct=prop.table(as.matrix(counts), margin=2)');
    ## Calculate the Jensen-Shannon Distance & parse output
    my $JSD_lines = $R->run('calc_JSD(ct)');
    my @JSD_split = split( /\n/, $JSD_lines );
    my $calc_JSD  = ( split /\s+/, $JSD_split[1] )[1];

    ## Calculate the JSdist CI for the model & parse output
    my $jsd_ci_lower;
    my $jsd_ci_upper;
    ## Bootstrap the model population while keeping the reference population_freq intact
    $R->run(
        'calc_JSD_boot_fxn = function (x_df, index) {
                            tmp_df <- data.frame( x_df[,1], x_df[index,2] )
                            return ( calc_JSD(tmp_df) ) 
                            }'
    );
    $R->run("JSD_boot <- boot(ct, calc_JSD_boot_fxn, R=$Bootstrap, parallel=\"multicore\", ncpus=$threads)");
    my $JSdist_ci = $R->run('boot.ci(JSD_boot, type="norm")');
    my $ci_data_line = ( split /\n/, $JSdist_ci )[8];
    if ( defined $ci_data_line ) {
        $ci_data_line =~ /\s+\((.+)\,\s+(.+)\)/;
        $jsd_ci_lower = $1;
        $jsd_ci_upper = $2;
    }
    else {
        $jsd_ci_lower = "NULL";
        $jsd_ci_upper = "NULL";
    }

    # Close R instance
    $R->stop();
    return $calc_JSD, $jsd_ci_lower, $jsd_ci_upper;
}

sub downsample_bam {
    my $desired_number_of_reads = shift;
    my $dir                     = shift;
    my $big_bam                 = shift;
    if ( !$desired_number_of_reads ) { confess "Error: Must pass &downsample_bam a desired_number_of_reads to downsample to.\n"; }
    if ( !$dir )                     { confess "Error: Must pass &downsample_bam a directory.\n"; }
    if ( !-e $big_bam )              { confess "Error: This bam doesn't exist: $big_bam\n"; }
    my $total  = wc($big_bam);
    my $prob   = $desired_number_of_reads / $total;
    my $header = run_cmd("samtools view -H $big_bam | head -n 1");
    if ( $header !~ /coordinate/ ) {
        run_cmd("samtools sort -@ $threads -m $sort_mem $big_bam $dir/downsample_psrt");
        $big_bam = "$dir/downsample_psrt.bam";
    }
GENERATE_BAM: run_cmd("$Picard DownsampleSam I=$big_bam O=$dir/downsampled.bam R=null P=$prob VALIDATION_STRINGENCY=LENIENT 2>>$dir/stderr.log");
    goto GENERATE_BAM if ( !-e "$dir/downsampled.bam" );
    my $downsampled_bam_count = wc("$dir/downsampled.bam");
    goto GENERATE_BAM if ( $downsampled_bam_count == 0 or $downsampled_bam_count % 2 );
    run_cmd("samtools sort -n -@ $threads -m 1G $dir/downsampled.bam $dir/nsort_downsampled 2>>$dir/stderr.log");
    return "$dir/nsort_downsampled.bam";
}

sub help {
    die "\nThis script will titrate the JSD for the [--reads] by downsampling a bam.
    --bam=              Bam to downsample from. *Mandatory*
    --pic=              Picard file for bam.
    --reads=            Range of reads to downsample to for calculating the JSD. ex: 2-10. Incremented by 2 over range.       
    --iterate=          [100]  Number of times iterate foreach [--reads]
    --boot=             [100]  Number of Bootstraps / iteration
    --threads=          [3]    Number of threads to sort and Bootstrap with. 
    --output_dir=       /path/for/output
    --output_prefix=    [--output_prefix]_jsd_titration.txt
    --sort_mem=         [500M] Amount of RAM to use / cpu to position sort (only used it unsorted bam passed).
    --sub_mem=          [7G]
    --sub_mail=         [0]
    --sub_name=         [JSDcalc]
    --help\n";
}

