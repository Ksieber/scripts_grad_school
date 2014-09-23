#!/usr/local/bin/perl
use warnings;
no warnings 'uninitialized';
use strict;
use run_cmd;
use Carp;
$Carp::MaxArgLen = 0;
use File::Basename;
use lib qw(/local/projects-t3/HLGT/scripts/lgtseek/lib/);
use LGTSeek;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,         'input|i=s',         'sort|S=i',       'window_size|w=i', 'step_size|s=i', 'min_cov|c=i',
    'min_quality|q=i', 'output_prefix|p=s', 'output_dir|o=s', 'tdo=i',           'Qsub|Q=i',      'region|r=s',
    'A=i',             'd=i',               'M_M=i',          'M_UM=i',          'help|?',
) or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} )        { help(); }
if ( !$options{input} )      { die "ERROR: Must give an input file with: --input= BAM or mpileup file." }
if ( !-e "$options{input}" ) { die "ERROR: Input file can't be found."; }
if ( !$options{window_size} || !$options{step_size} ) { die "ERROR: Must pass the --window_size & --step_size. Use --Help for more info.\n"; }

my $lgtseq = LGTSeek->new2( \%options );
our $window_size = defined $options{window_size} ? $options{window_size}   : "10";
our $step_size   = defined $options{step_size}   ? $options{step_size} - 1 : "5";
our $min_cov     = defined $options{min_cov}     ? "$options{min_cov}"     : "1";
our $tdo         = defined $options{tdo}         ? "$options{tdo}"         : "0";

my $sort = defined $options{sort} ? "$options{sort}" : "0";
my $qsub = defined $options{Qsub} ? "$options{Qsub}" : "0";

my $region  = defined $options{region}      ? "'$options{region}'"       : undef;
my $d       = defined $options{d}           ? $options{d}                : 100000;
my $quality = defined $options{min_quality} ? "-q $options{min_quality}" : "-q 0";
my $view    = "$quality -hu";                              ## Default
my $A       = defined $options{A} ? "$options{A}" : "1";
my $mpileup = ( $A == 0 ) ? "-d $d" : "-Ad $d";            ## Default

my $window;                                                # keys: coverage (total coverage seen in window), start (actual start), stop (expected stop), bases (total in window seen), chr
my $next_window;
my $fh = open_input( $options{input} );

while (<$fh>) {
    chomp;
    my $line = $_;
    my ( $chr, $pos, $coverage ) = ( split /\t/, $line )[ 0, 1, 3 ];
    if ( $chr ne $window->{'chr'} ) {
        if ( defined $window->{'chr'} ) { print_coverage($window); }
        $window = init_window($line);
        next;
    }

    # Initialize the next window
    if ( $pos == ( $window->{'start'} + $step_size ) ) { $next_window = init_window($line); }

    # Extend the NEXT window
    if ( $pos > $next_window->{'start'} && $pos < $next_window->{'stop'} ) {
        $next_window->{'coverage'} += $coverage;
        $next_window->{'bases'}++;
    }

    # Extend the CURRENT window
    if ( $pos > $window->{'start'} && $pos < $window->{'stop'} ) {
        $window->{'coverage'} += $coverage;
        $window->{'bases'}++;
        next;
    }

    # If past the current window size, print current window and move to next window.
    if ( $pos >= $window->{'stop'} ) {
        print_coverage($window);
        if ( $pos > $next_window->{'stop'} ) { $next_window = undef; }
        $window      = $next_window;
        $next_window = init_window($line);
    }
}

print_coverage($window);
close $fh;

sub open_input {
    my $raw_input = shift;
    my ( $fn, $path, $suffix ) = fileparse( $raw_input, ( @{ $lgtseq->{mpileup_suffix_list} }, @{ $lgtseq->{bam_suffix_list} } ) );

    if ( $lgtseq->empty_chk( { input => $raw_input } ) == 1 ) { confess "ERROR: Input file:$raw_input is empty."; }
    ## Filehandle to return.
    my $fh;
    ## Mpileup input
    map {
        if ( $suffix eq $_ ) {
            ## Warn we can't filter quality with mpileup input.
            if ( defined $options{min_quality} ) {
                print STDERR
                    "\n****** Warning ****** \nInput: $options{input} is a mpileup text file. --min_quality needs a bam file input to filter for min_quality. Proceeding with no min_quality filtering.\n*** Warning *** \n\n";
            }
            ## Open input file.
            open( $fh, "<", "$raw_input" ) || confess "ERROR: Can't open input: $raw_input because: $!\n";
            return $fh;
        }
    } @{ $lgtseq->{mpileup_suffix_list} };

    ## Position sorted bam
    my @psort_bams = ( '_pos-sort.bam', '_psort.bam', 'srt.bam' );
    map {
        if ( $suffix eq $_ ) {
            print STDERR "FOOBAR\n";
            open( $fh, "-|", "samtools view $view $raw_input $region | samtools mpileup $mpileup - " )
                || confess "ERROR: Can't open input: $raw_input because: $!\n";
            return $fh;
        }
    } @psort_bams;

    ## Name sorted bam || --sort=1
    map {
        if ( $suffix eq $_ || $sort == 1 ) {
            print STDERR "Sorting input bam . . .\n";
            open( $fh, "-|", "samtools sort -m 1G -o $raw_input - | samtools view $view - $region | samtools mpileup $mpileup - " )
                || confess "ERROR: Can't open input: $raw_input because: $!\n";
            return $fh;
        }
    } @{ $lgtseq->{bam_suffix_list} };
}

sub print_coverage {
    my $window          = shift;
    my $window_coverage = sprintf( "%.2f", $window->{'coverage'} / $window->{'bases'} );
    my $actual_stop     = ( $window->{'start'} + $window->{'bases'} );
    my $step_pos        = $window->{'start'} + $step_size;
    my $to_print        = $tdo == 1 ? "$window_coverage\t$window->{chr}\t$step_pos\n" : "$window_coverage\t$window->{chr}\:$window->{start}\-$actual_stop\n";
    print $to_print if ( $window_coverage >= $min_cov );
}

sub init_window {
    my $line = shift;
    my ( $chr, $pos, $coverage ) = ( split( /\s/, $line ) )[ 0, 1, 3 ];
    my %window;
    $window{chr}      = $chr;
    $window{start}    = $pos;
    $window{stop}     = $pos + $window_size;
    $window{coverage} = $coverage;
    $window{bases}    = 1;
    return \%window;
}

sub help {
    die "Help: This script will calculate coverage with a sliding window. 
    --input|i=             Mpileup or bam.
        --sort|S=          <0|1> [0] 1= Position sort input bam.
        --window_size|w=   Size of the window to calculate coverage for. 
        --step_size|s=     Size of the steps inbetween windows.
        --min_cov|c=       Minimal coverage to report.
        --min_quality|q=   [0] Minimal quality score for mapping.
        --output_prefix=   Prefix for the output.
        --output_dir|o=    Directory for output.
        --tdo=             <0|1> [0] 1= Tab Delimited Output of chr and position. 0= IGV style.
        --Qsub|Q=          <0|1> [0] 1= Qsub job to grid.
    --help|?\n";
}

