#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
my $VERSION = "2.01";
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

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,   'pic1=s',     'pic2=s', 'bam1=s',      'bam2=s',   'output_dir|o=s', 'output_prefix|p=s', 'ref=s',
    'subdirs=i', 'sort_mem=s', 'boot=i', 'threads|t=i', 'Qsub|q=i', 'help|?',         'sub_name=s',        'project=s',
) or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} ) { &help; }
if ( !$options{pic1} and !$options{bam1} ) { die "Error: Must pass input with --pic1 or --bam1. Please try again.\n"; }
if ( !$options{pic2} and !$options{bam2} ) { die "Error: Must pass input with --pic1 or --bam1. Please try again.\n"; }

my $lgtseq    = LGTSeek->new2();
my $Picard    = "$lgtseq->{java_bin} \-$lgtseq->{java_opts} -jar $lgtseq->{Picard_jar}";
my $threads   = defined $options{threads} ? $options{threads} : "3";
my $sort_mem  = defined $options{sort_mem} ? $options{sort_mem} : "1G";
my $Bootstrap = defined $options{boot} ? $options{boot} : "1000";
my ( $fn, $path, $suf )
    = ( defined $options{pic1} and -e $options{pic1} )
    ? fileparse( $options{pic1}, $lgtseq->{suffix_regex} )
    : fileparse( $options{bam1}, $lgtseq->{bam_suffix_list} );
my $output_dir = defined $options{output_dir} ? $options{output_dir} : "$path/tmp";
mk_dir("$output_dir");
my $output_prefix = defined $options{output_prefix} ? $options{output_prefix} : $fn;

my $picard1 = ( defined $options{pic1} and -e $options{pic1} ) ? $options{pic1} : &bam2picard( $options{bam1} );    # Returns picard_file path
my $picard2 = ( defined $options{pic2} and -e $options{pic2} ) ? $options{pic2} : &bam2picard( $options{bam2} );

my $pop1 = &pop_insert_size_from_picard_file($picard1);                                                             # Returns Hash_ref->{insert_size}=count
my $pop2 = &pop_insert_size_from_picard_file($picard2);

my ( $jsd, $jsd_ci_min, $jsd_ci_max ) = &calc_jsd( $pop1, $pop2 );

printf STDOUT ( "%-20s%-20s%-20s\n",       "JSD", "CI_min",    "CI_max" );
printf STDOUT ( "%-20.4f%-20.4f%-20.4f\n", $jsd,  $jsd_ci_min, $jsd_ci_max );

# Clean up tmp dir if no [output_dir] was passed
if ( !$options{output_dir} and -e "$path/tmp/" ) { run_cmd("rm -rf $output_dir"); }

sub bam2picard {
    my $bam = shift;

    my $header = run_cmd("samtools view -H $bam | head -n 1");
    if ( $header !~ /coordinate/ ) {
        run_cmd("samtools sort -@ $threads -m $sort_mem $bam $output_dir/$output_prefix\.srt");
        $bam = "$output_dir/$output_prefix\.srt.bam";
    }
    run_cmd("$Picard CollectInsertSizeMetrics AS=true I=$bam H=/dev/null O=$output_dir/$output_prefix\_insert.metrics M=0 VALIDATION_STRINGENCY=SILENT 2>>$output_dir/picard_stderr.log");
    return "$output_dir/$output_prefix\_insert.metrics";
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
    $R->run('pop_1_count   = numeric()');
    $R->run('pop_2_count = numeric()');
    if ( scalar( keys %$pop1 ) > scalar( keys %$pop2 ) ) {
        foreach my $insert ( sort { $a <=> $b } keys %$pop1 ) {
            my $pop_1_isize_count = $pop1->{$insert};
            my $pop_2_isize_count = defined $pop2->{$insert} ? $pop2->{$insert} : "0.000000000001";
            ## Add the pop and model count of each insert size to the matrix by col.
            $R->run("pop_1_count   = c( pop_1_count,   $pop_1_isize_count   )");
            $R->run("pop_2_count = c( pop_2_count, $pop_2_isize_count )");
        }
    }
    else {
        foreach my $insert ( sort { $a <=> $b } keys %$pop2 ) {
            my $pop_1_isize_count = defined $pop1->{$insert} ? $pop1->{$insert} : "0.000000000001"; 
            my $pop_2_isize_count = defined $pop2->{$insert} ? $pop2->{$insert} : "0.000000000001";
            ## Add the pop and model count of each insert size to the matrix by col.
            $R->run("pop_1_count   = c( pop_1_count,   $pop_1_isize_count   )");
            $R->run("pop_2_count = c( pop_2_count, $pop_2_isize_count )");
        }
    }

    $R->run('counts=data.frame(pop_1_count,pop_2_count)');
    ## Calculate the proportion of each Isize in respective populations(ref & model)
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

sub help {
    die "This script will calculate the JSD between to populations of reads.
    --pic1=             Picard file for population 1.          (Faster than bam input)
    --pic2=             Picard file for population 2.          (Faster than bam input)
    --bam1=             Bam for population1. (Much faster if already position sorted.)
    --bam2=             Bam for population2. (Much faster if already position sorted.)
    --threads=          [3] < # > Number of threads to sort and Bootstrap with. 
    --sort_mem=         [1G] Amount of RAM to use / cpu to position sort (only used it unsorted bam passed).
    --boot=             [1000]  Number of Bootstraps
    --output_dir=   
    --output_prefix=
    --help\n";
}
