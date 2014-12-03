#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
my $VERSION = "2.00";
use warnings;
use strict;
use Carp;
use run_cmd;
use print_call;
use mk_dir;
use bwa;
use read_bam;
use empty_chk;
use File::Basename;
use Bio::DB::Fasta;
use POSIX;
use Statistics::R;
use Math::NumberCruncher;
use Data::Dumper;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,     'lgt_ref|i=s', 'input|i=s',   'read_length=i', 'iterate=i', 'reads=s',      'output_dir|o=s', 'output_prefix|p=s', 'ref=s', 'subdirs=i',
    'threads|t=i', 'Qsub|q=i',    'wgsim_bin=s', 'help|?',        'chr1=s',    'ref1_start=i', 'sub_name=s',
) or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} ) { &help; }

# Check to make sure we have all the options we will need
if ( !$options{input} ) { die "Error: Must pass a --input=<Picard File> Please try again.\n"; }
if ( !$options{reads} ) { die "Error: Must pass a --reads= number of random reads to generate to model JSD. Please try again\n"; }
if ( !$options{ref} and !$options{lgt_ref} ) { die "Error: Must pass a ref with --ref=standard ref. or --lgt_ref=merge2lgtbam.fa. Please try again.\n"; }

# Sumbit the job to the grid
if ( defined $options{Qsub} and $options{Qsub} == 1 ) {
    if ( !$options{sub_name} ) { $options{sub_name} = "jsdCalc"; }
    Qsub_script( \%options );
}

# Set default and global variables
my @reads_number_to_process = -e $options{reads}            ? @{ &read_in_list( $options{reads} ) } : split( /,/, join( ',', $options{reads} ) );
my $wgsim                   = defined $options{wgsim_bin}   ? $options{wgsim_bin}                   : "/home/ksieber/bin/wgsim";
my $read_length             = defined $options{read_length} ? $options{read_length}                 : "51";
my $iterate                 = defined $options{iterate}     ? $options{iterate}                     : "1000";
my $threads                 = defined $options{threads}     ? $options{threads}                     : "8";
my $ref1                    = $options{ref};
my ( $fn, $path, $suf ) = fileparse( $options{input}, ( ".std_insert.metrics", ".metrics", ".txt", qr/\.[^\.]+/ ) );
my $overall_output_dir = defined $options{output_dir}    ? $options{output_dir}    : $path;
my $output_prefix      = defined $options{output_prefix} ? $options{output_prefix} : $fn;
if ( defined $options{subdirs} and $options{subdirs} == 1 ) { $overall_output_dir = "$overall_output_dir/$output_prefix/"; }
mk_dir("$overall_output_dir");
print_notebook( \%options );

# First, get the insert size and deviation of the library from the Picard file.
my $insert_size;
my $stdev;
my $picard_header = run_cmd("head -n 8 $options{input}");
my @picard_lines = split( /\n/, $picard_header );
if   ( $picard_lines[6] !~ /^MEDIAN_INSERT_SIZE/ ) { confess "Error: The [--input] doesn't look right. Please fix it and try agian.\n"; }
else                                               { ( $insert_size, $stdev ) = ( split /\t/, $picard_lines[7] )[ 0, 1 ]; }

my @jsd_values;
my @jsd_ci_min_values;
my @jsd_ci_max_values;
my $n = 1;

my $generated_ref = &generate_standard_ref( $overall_output_dir, $ref1 ) if ( $options{ref} and !$options{lgt_ref} );
my ( $lgt_ref1, $lgt_ref2 ) = &generate_lgt_refs( $overall_output_dir, $options{lgt_ref}, $ref1, $options{ref2} ) if ( defined $options{lgt_ref} and -e $options{lgt_ref} );

foreach my $number_of_reads (@reads_number_to_process) {
    my $subdir_name = "$number_of_reads\_reads";
    my $output_dir = $overall_output_dir =~ /$number_of_reads\_reads$/ ? $overall_output_dir : "$overall_output_dir/$subdir_name/";
    mk_dir($output_dir);
    open( my $OUT, ">", "$output_dir/JSD_values_conf-intrvls.txt" );
    print $OUT "\# titrate_jsd_random_reads_from_INT.pl\tVersion:$VERSION\n";
    print $OUT "\# iterate=$iterate\treads=" . join( ",", @reads_number_to_process ) . "\n";
    print $OUT "\# Picard_input=$options{input}";
    if ( defined $options{ref} )     { print $OUT "\tref1=$ref1"; }
    if ( defined $options{lgt_ref} ) { print $OUT "\tlgt_ref=$options{lgt_ref}"; }
    print $OUT "\n";

    # Second, start iterating for [--iterate] number of times for each $number_of_reads
    for ( $n = 1; $n <= $iterate; $n++ ) {
        &print_progress($n);
        printf $OUT ( "%-20.5f", $n );
        my $tmp_dir = "$output_dir/tmp_$fn";
        mk_dir("$tmp_dir");

        my $wgsim_bam
            = ( defined $generated_ref and -e $generated_ref )
            ? &process_standard_reference( $number_of_reads, $tmp_dir, $generated_ref )
            : &process_merge2lgt_reference( $number_of_reads, $tmp_dir, $lgt_ref1, $lgt_ref2 );

        #  Prepare data for and calculate JSD
        ## Build model population count
        my %model;
        my $BAM = open_bam($wgsim_bam);
        while ( my $read1 = read_bam($BAM) ) {
            my $read2 = read_bam($BAM);
            $model{ abs( $read2->{insert} ) }++;
        }
        ## Build reference population count
        my %reference;
        open( PIC, "<", "$options{input}" ) or confess "Error: Unable to open the [--input]\n";
        my @header              = <PIC> . <PIC> . <PIC> . <PIC> . <PIC> . <PIC> . <PIC>;
        my @picard_insert_sizes = <PIC> . <PIC> . <PIC>;
        my $foobar              = <PIC> . <PIC> . <PIC>;
        while (<PIC>) {
            next if ( $_ !~ /^\d+/ );
            my ( $i_size, $fr, $rf, $tandem ) = split( /\t/, $_ );
            $reference{$i_size} = $fr;
        }
        close PIC;
        ## Calculate JSD
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
        foreach my $key ( sort { $a <=> $b } keys %reference ) {
            my $pop_isize_count = $reference{$key};
            my $model_isize_count = defined $model{$key} ? $model{$key} : "0.0000001";
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
        push( @jsd_values, $calc_JSD );
        printf $OUT ( "%-20.5f", $calc_JSD );
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
        $R->run("JSD_boot <- boot(ct, calc_JSD_boot_fxn, R=100, parallel=\"multicore\", ncpus=$threads)");
        my $JSdist_ci = $R->run('boot.ci(JSD_boot, type="norm")');
        my $ci_data_line = ( split /\n/, $JSdist_ci )[8];
        if ( defined $ci_data_line ) {
            $ci_data_line =~ /\s+\((.+)\,\s+(.+)\)/;
            $jsd_ci_lower = $1;
            $jsd_ci_upper = $2;
            push( @jsd_ci_min_values, $jsd_ci_lower );
            push( @jsd_ci_max_values, $jsd_ci_upper );
            printf $OUT ( "%-20.4f%-20.4f", $jsd_ci_lower, $jsd_ci_upper );    ##
        }
        else {
            $jsd_ci_lower = "NULL";
            $jsd_ci_upper = "NULL";
            push( @jsd_ci_min_values, $jsd_ci_lower );
            push( @jsd_ci_max_values, $jsd_ci_upper );
            printf $OUT ( "%-20s%-20s", $jsd_ci_lower, $jsd_ci_upper );        ##
        }
        printf $OUT "\n";

        # Close R instance
        $R->stop();
        run_cmd("rm -rf $tmp_dir");
    }

    # After [--iterate] is complete for the given $number_of_reads, summarize the $number_of_reads data.
    close $OUT;
    open( my $TXT, ">", "$output_dir/JSD_titration_summary.txt" ) or confess "Error: Unable to open the JSD summary.txt file: $output_dir/JSD_avg_ranges.txt\n";
    &print_summary_stats( $TXT, "JSD_values", \@jsd_values );
    &print_summary_stats( $TXT, "JSD_ci-min", \@jsd_ci_min_values );
    &print_summary_stats( $TXT, "JSD_ci-max", \@jsd_ci_max_values );
    close $TXT;
}
if ( defined $options{lgt_ref} and -e $options{lgt_ref} ) {
    run_cmd("rm $lgt_ref1 $lgt_ref2");
}
else {
    run_cmd("rm $generated_ref");
}

print_complete( \%options );

## COMPLETE

## Subroutines below:
sub print_summary_stats {
    my $FH        = shift;
    my $line_name = shift;
    my $array     = shift;
    my ( $max, $min ) = Math::NumberCruncher::Range($array);
    my $avg       = Math::NumberCruncher::Mean($array);
    my $deviation = Math::NumberCruncher::StandardDeviation($array);
    printf $FH ( "%-20s", $line_name );
    printf $FH ( "%-20.6f%-20.6f%-20.6f%-20.6f", $avg, $deviation, $min, $max );
    print $FH "\n";
}

sub process_standard_reference {
    my $number_of_reads = shift;
    my $tmp_dir         = shift;
    my $generated_ref   = shift;

    #  Generate random PE reads from the reference
    run_cmd("$wgsim -e 0 -d $stdev -s $insert_size -N $number_of_reads -1 $read_length -2 $read_length $generated_ref $tmp_dir/Test_1.fq $tmp_dir/Test_2.fq 2>>$tmp_dir/wgsim_stderr.txt");
    if ( empty_chk("$tmp_dir/Test_1.fq") == 1 or empty_chk("$tmp_dir/Test_2.fq") == 1 ) {
        confess "Error: The reference isn't large enough to generate wgsim random reads with the [--input] insert size and stdev.\n";
    }
    ## Map the random reads to the ref1
    my $wgsim_bam = bwa_align(
        "$tmp_dir/Test_1.fq,$tmp_dir/Test_2.fq",
        $ref1,
        {   output_dir    => $tmp_dir,
            output_prefix => "test_$n",
            threads       => $threads,
        }
    );
    ## Check to make sure we had reads map to the reference
    if ( empty_chk($wgsim_bam) == 1 ) {
        print STDERR "Error: No reads mapped to the random_ref.\n";
        run_cmd("rm -rf $tmp_dir");
        return undef;
    }
    return $wgsim_bam;
}

sub process_merge2lgt_reference {
    my $number_of_reads = shift;
    my $tmp_dir         = shift;
    my $lgt_ref1        = shift;
    my $lgt_ref2        = shift;

    #  Generate random SE reads from the lgt_Ref1
    run_cmd("$wgsim -e 0 -d 0 -s $stdev -N $number_of_reads -1 $read_length $lgt_ref1 $tmp_dir/Test_1.fq /dev/null 2>>$tmp_dir/wgsim_stderr.txt");
    if ( empty_chk("$tmp_dir/Test_1.fq") == 1 ) {
        confess "Error: Unable to generate SE reads for $lgt_ref1.\nThe reference is probably too short for the STDEV and read length. Check $tmp_dir/wgsim_stderr.txt\n";
    }
    ##  Generate random SE reads from the lgt_Ref2
    run_cmd("$wgsim -e 0 -d 0 -s $stdev -N $number_of_reads -2 $read_length $lgt_ref2 /dev/null $tmp_dir/Test_2.fq 2>>$tmp_dir/wgsim_stderr.txt");
    if ( empty_chk("$tmp_dir/Test_2.fq") == 1 ) {
        confess "Error: Unable to generate SE reads for $lgt_ref2.\nThe reference is probably too short for the STDEV and read length. Check $tmp_dir/wgsim_stderr.txt\n";
    }
    ## Map the random reads to the ref1
    my $wgsim_bam = bwa_align(
        "$tmp_dir/Test_1.fq,$tmp_dir/Test_2.fq",
        $options{lgt_ref},
        {   output_dir    => $tmp_dir,
            output_prefix => "test_$n",
            threads       => $threads,
        }
    );
    ## Check to make sure we had reads map to the reference
    if ( empty_chk($wgsim_bam) == 1 ) {
        print STDERR "Error: No reads mapped to the random_ref.\n";
        run_cmd("rm -rf $tmp_dir");
        return undef;
    }
    return $wgsim_bam;
}

sub generate_standard_ref {
    my $output_dir = shift;
    my $reference  = shift;
    if ( !-e $reference ) { confess "Error: Reference doesn't exist: $reference\n"; }

    # Calculate the min. sequence length needed
    my $ref_size = ( $insert_size * 3 ) + $stdev + 2;

    # Pull the reference sequence
    my $ref_db = Bio::DB::Fasta->new($ref1);
    my $chr = defined $options{chr1} ? $options{chr1} : run_cmd("head -n 1 $ref1 | cut -f1 -d \' \'");
    $chr =~ s/^\>//;
    my $start = defined $options{ref1_start} ? $options{ref1_start} : "0";
    my $end   = $start + $ref_size;
    my $seq   = $ref_db->seq( $chr, $start => $end );

    # Print the selected reference sequence.
    open( my $FA, ">", "$output_dir/Generated.fa" ) or confess "Error: Unable to open the Generated.fa: $output_dir/Generated.fa\n";
    print $FA "\>$chr\n";
    print $FA "$seq";
    close $FA;

    # Complete
    return "$output_dir/Generated.fa";
}

sub generate_lgt_refs {
    my $output_dir     = shift;
    my $merged2lgt_ref = shift;
    my $ref1           = shift;
    my $ref2           = shift;
    my ( $merged2lgt_fn, $merged2lgt_path, $merged2lgt_suff ) = fileparse( $merged2lgt_ref, ( ".fa", ".fasta" ) );
    if ( -e "$merged2lgt_path/Ref_img_cords.txt" ) {
        open( my $IN, "<", "$merged2lgt_path/Ref_img_cords.txt" ) or confess "Error: Unable to open: $merged2lgt_path/Ref_img_cords.txt\n";
        ## Moving past the header
        my $null_header = <$IN>;
        ## lgt ref1 cords
        my $lgt_ref1_cords = ( split /\t/, <$IN> )[1];
        my ( $start1, $end1 ) = ( split /-/, $lgt_ref1_cords )[ 0, 1 ];
        ## N string, not needed
        my $null = <$IN>;
        ## lgt ref2 cords
        my $lgt_ref2_cords = ( split /\t/, <$IN> )[1];
        my ( $start2, $end2 ) = ( split /-/, $lgt_ref2_cords )[ 0, 1 ];
        close $IN;
        ## Grab lgt1 & 2 sequences
        my $chr           = run_cmd("head -n1 $merged2lgt_ref | cut -f2 -d '>'");
        my $merged2lgt_db = Bio::DB::Fasta->new($merged2lgt_ref);
        my $seq1          = $merged2lgt_db->seq( $chr, $start1 => $end1 );
        my $seq2          = $merged2lgt_db->seq( $chr, $start2 => $end2 );
        ##  Build lgt_ref1 fasta ref.
        open( my $FA1, ">", "$output_dir/lgt_ref1.fa" ) or confess "Error: Unable to open output: $output_dir/lgt_ref1.fa\n";
        print $FA1 ">lgt_ref1\n";
        print $FA1 "$seq1";
        close $FA1;
        ##  Build lgt_ref2 fasta ref.
        open( my $FA2, ">", "$output_dir/lgt_ref2.fa" ) or confess "Error: Unable to open output: $output_dir/lgt_ref2.fa\n";
        print $FA2 ">lgt_ref2\n";
        print $FA2 "$seq2";
        close $FA2;
        return "$output_dir/lgt_ref1.fa", "$output_dir/lgt_ref2.fa";
    }
    else {
        confess "derp derp";
    }
}

sub print_progress {
    my $n = shift;
    print STDERR "Progress . . . \n" if ( $n eq 1 );
    print STDERR "  ."               if ( $n =~ /\d*25$/ );
    print STDERR "  ."               if ( $n =~ /\d*50$/ );
    print STDERR "  ."               if ( $n =~ /\d*75$/ );
    print STDERR "  $n\n"            if ( $n eq 100 );
    print STDERR "  $n\n"            if ( $n eq 200 );
    print STDERR "  $n\n"            if ( $n eq 300 );
    print STDERR "  $n\n"            if ( $n eq 400 );
    print STDERR "  $n\n"            if ( $n eq 500 );
    print STDERR "  $n\n"            if ( $n eq 600 );
    print STDERR "  $n\n"            if ( $n eq 700 );
    print STDERR "  $n\n"            if ( $n eq 800 );
    print STDERR "  $n\n"            if ( $n eq 900 );
}

sub help {
    die "   This script will take a INT reference and generate [--reads=] number of reads with the given library insert size and distribution.
    Using these random reads, the optimal Jenson-Shannon Distance will be calculated.
    This process is repeated [--iterate] number of times to access the best possible JSD for [--reads] number of reads.
        --input|P=              Picard insert metrics for the original and entire library
        --reads=                The number of reads to generate for calculating the JSD. A single #, comma seperated string, or a file w/ 1 # /line.
        --ref=                  Standard Fasta Reference. 
        --lgt_ref=              Fasta reference from merge_2lgt_bams.pl ## Not implemented yet.
        --iterate=                  Number of times to iterate over each number of reads.
        --output_dir=               /path/for/output/
        --threads=   
         \n";
}
