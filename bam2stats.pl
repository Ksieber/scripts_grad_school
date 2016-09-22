#!/usr/bin/perl
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/', '/local/projects-t3/HLGT/scripts/lgtseek/lib/' );
use strict;
use warnings;
use Time::SoFar;
use File::Basename;
use LGTSeek;
use run_cmd;
use mk_dir;
use setup_input;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
our $results = GetOptions( \%options, 'input|i=s', 'input_list|I=s', 'output_dir|o=s', 'subdirs=s', 'append=s', 'Qsub|q=s', 'sub_mem=s', 'sub_wd=s', 'help|?' )
    or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} ) {
    die "This script will calculate the number of reads in a bam: Total, MM, MU, UU, SC
    --input|i=              <BAM> (Also accepts single bam from \$ARGV[0])
    --input_list|I=         <list of BAMS>
    --append=               </path/file/append.txt> Add output stats to this file.
    --output_dir|o=         </path/for/output/>
    --subdirs=              <0|1> [0] 1= Put output into a subdirectory under --output_dir.
    --Qsub|q=               <0|1> [0] 1= Qsub the script.
    --sub_wd=               [--output_dir] /dir/for/qsub.e# files 
    --help|?\n";
}

if ( !$options{input} && $ARGV[0] ) { $options{input} = $ARGV[0]; }
if ( !$options{input} && !$options{input_list} ) { die "Must give an --input=<BAM> or --input_list=<list_of_bams>\n"; }
my $append = defined $options{append} ? "$options{append}" : "0";

my $lgtseek = LGTSeek->new2( \%options );

my $input = setup_input( \%options );

foreach my $input (@$input) {
    if ( $lgtseek->empty_chk( { input => $input } ) == 1 ) { print STDERR "This is an empty bam! $input\n"; next; }
    if ( $lgtseek->{Qsub} == 1 ) {
        ## If we are in the orignal call, change input from list to a single file
        if ( $options{input_list} ) { $options{input} = $input; }
        ## Build qsub command
        my $cmd = "/home/ksieber/scripts/bam2stats.pl";
        foreach my $key ( keys %options ) {
            next if ( $key =~ /Qsub/ );
            if ( defined $options{$key} ) { $cmd = $cmd . " --$key=$options{$key}" }
        }
        my $output_dir;
        if ( defined $options{append} ) {
            my ( $file, $dir, $foo ) = fileparse( $options{append} );
            $output_dir = $dir;
        }
        elsif ( defined $options{output_dir} ) {
            $output_dir = $options{output_dir};
        }
        else {
            $output_dir = run_cmd("pwd");
        }
        mk_dir($output_dir);
        my $wd_dir = defined $options{sub_wd} ? $options{sub_wd} : $output_dir;
        mk_dir($wd_dir);

        ## submit command to grid
        Qsub(
            {   cmd      => "$cmd",
                sub_name => "bam2stats",
                sub_mem  => "1G",
                wd       => $wd_dir,
            }
        );
        ## Skip to next input for qsub
        die "+++ Finished submiting job(s) to the grid. +++\n";
    }
    my $start = Time::SoFar::runinterval();
    open( my $in, "-|", "samtools view $input" ) or die "Can't open input: $input\n";
    ##
    my @read_types = ( 'Total', 'MM', 'MU', 'UU', 'SC' );
    my %counts;
    map { $counts{$_} = 0; } @read_types;
    while (<$in>) {
        chomp;
        my ( $flag, $cigar ) = ( split /\t/, $_ )[ 1, 5 ];
        my $converted_flag = $lgtseek->_parseFlag($flag);
        if ( $cigar =~ /(\d+)M(\d+)S/        && $2 >= 24 )                        { $counts{SC}++; }
        if ( $cigar =~ /(\d+)S(\d+)M/        && $1 >= 24 )                        { $counts{SC}++; }
        if ( !$converted_flag->{'qunmapped'} && !$converted_flag->{'munmapped'} ) { $counts{MM}++; }
        if ( !$converted_flag->{'qunmapped'} && $converted_flag->{'munmapped'} )  { $counts{MU}++; }
        if ( $converted_flag->{'qunmapped'}  && !$converted_flag->{'munmapped'} ) { $counts{MU}++; }
        if ( $converted_flag->{'qunmapped'}  && $converted_flag->{'munmapped'} )  { $counts{UU}++; }
        $counts{Total}++;
    }
    ##
    my $OUT;
    my $print_header;
    if ( $options{output_dir} ) {
        my ( $fn, $path, $suf ) = fileparse( $input, ".bam" );
        my $out = $options{output_dir} ? "$options{output_dir}/$fn.stats" : "$path$fn.stats";
        open( $OUT, "> $out" ) or die "Error: Can not open: $out\n";
        $print_header = 1;
    }
    elsif ( $options{append} ) {
        if   ( -e $options{append} ) { $print_header = 0; }
        else                         { $print_header = 1; }
        open( $OUT, ">> $options{append}" ) or die "Error: Can not open: $append\n";
    }
    else {
        $OUT          = *STDOUT;
        $print_header = 1;
    }
    ##
    if ( $print_header == 1 ) { print $OUT "BAM\t"; print $OUT join( "\t", @read_types ) . "\n"; }
    ##
    print $OUT "$input\t";
    map { print $OUT "$counts{$_}\t" } @read_types;
    my $finished = Time::SoFar::runinterval();
    print $OUT "Time elapsed: $finished\n";
}

