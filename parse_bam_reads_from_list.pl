#!/usr/local/bin/perl
use strict;
use warnings;
use run_cmd;
use mk_dir;
use File::Basename;
use Getopt::Long qw (:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'good_list=s', 'bad_list=s', 'output_prefix|p=s', 'output_dir|o=s', 'sub_name=s', 'Qsub|q=i', 'help|?' )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) {
    die "HELP: This script will parse an bam file for the desired reads.
--input|i=            A bam file to pull the reads from.
--good_list=          A list of reads that you want from the bam file.
--bad_list=           A list of reads to remove from the bam file.
--output_prefix=     \$prefix.bam
--output_dir|o=       Directory for output
--Qsub|q=             <0|1> [0] 1= Submit the job to the grid.\n";
}

if ( !$options{input} ) {
    die "ERROR: You must give a bam file to parse. Use --input=<FILE>\n";
}
if ( !$options{good_list} && !$options{bad_list} ) {
    die "ERROR: Must give a list reads to pull from the .bam file. Use either: --good_list=<FILE> or --bad_list=<FILE>\n";
}

if ( defined $options{Qsub} and $options{Qsub} == 1 ) { $options{sub_name} = "parsebam"; Qsub_script( \%options ); }
my ( $file, $path, $suf ) = fileparse( $options{input}, ".bam" );
my $prefix = $options{output_prefix} ? $options{output_prefix} : "$file\_filtered";
my $dir    = $options{output_dir}    ? $options{output_dir}    : $path;
mk_dir($dir);
my $out = "$dir\/$prefix";

my %mapped_reads;    ## This is a terrible name for the hash. It is used to keep track of the reads to parse on. Ie good or bad reads depending on what is given.

if ( $options{good_list} || $options{bad_list} ) {
    my $in;
    if ( $options{good_list} ) {
        print STDERR "Opening GOOD list: $options{good_list} to read in the desired reads . . .\n";
        my $wc = `wc -l $options{good_list}`;
        my @wc = split( /\s/, $wc );
        if ( $wc[0] eq 0 ) {
            die "ERROR: The good list: $options{good_list} is empty! Cannot pull reads from this list.\n";
        }
        open( IN, "<", $options{good_list} ) or die "Unable to open input good_list: $options{good_list}\n";
    }
    if ( $options{bad_list} ) {
        print STDERR "Opening BAD list: $options{bad_list} to read in the reads to remove . . .\n";
        my $wc = `wc -l $options{bad_list}`;
        my @wc = split( /\s/, $wc );
        if ( $wc[0] eq 0 ) {
            die "ERROR: The good list: $options{good_list} is empty! Cannot ignore these reads.\n";
        }
        open( IN, "<", $options{bad_list} ) or die "Unable to open input bad_list: $options{bad_list}\n";
    }
    while (<IN>) {
        chomp;

        #      $_=~s/_(1|2)$//;  This line will remove _1 or _2 if the reads have it. Maybe want to use this depending... USE WITH CAUTION.
        $mapped_reads{$_}++;
    }
    print STDERR "Finished reading in list of reads to parse.\n";
}
close IN;

open( OUT, "| samtools view -S - -bo $out\.bam" ) or die "Unable to open output: $out\.bam because: $!\n";
print STDERR "Opening output: $out\.bam\n";
open( BAM, "-|", "samtools view -h $options{input}" ) or die "Unable to run samtools view on $options{input} because: $!\n";
print STDERR "Starting to read through: $options{input}  parsing for the desired reads . . .\n";
while (<BAM>) {
    if ( $_ =~ /^@/ ) { print OUT; next; }    ## Print header lines
    my @f = split( /\t/, $_ );
    my $read = $f[0];

    #print STDERR "$read\n";      ## Working
    #    $read=~s/_(1|2)$//;  This line will remove _1 or _2 if the reads have it. Maybe want to use this depending... USE WITH CAUTION.
    if ( $options{good_list} && $mapped_reads{$read} ) {
        print OUT;
    }
    if ( $options{bad_list} ) {
        next if ( $mapped_reads{$read} );
        print OUT;
    }
}
print STDERR "Finished parsing the BAM file.\n";
close BAM;
close OUT;
print STDERR "Finished writing: $out\.bam\n";

