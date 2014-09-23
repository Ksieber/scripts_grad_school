#!/usr/bin/perl -I /home/ksieber/scripts/ -I /home/ksieber/perl5/lib/perl5/
use warnings;
no warnings 'uninitialized';
use strict;
use run_cmd;
use setup_input;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input=s', 'input_list=s', 'sort=s', 'output_prefix=s', 'output_dir=s', 'Qsub=s', 'region=s', 'A=s', 'd=s', 'M_M=s', 'M_UM=s', 'help', 'ref=s', )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) {
    die "HELP: This script will take a bam file and calculate the coverage.
        --input=            Input bam (Position sorted already)
        --input_list=       List of bams to process.
        --sort=             <0|1> [0] 1= Sort input by position. In order to calc. mpileup input must be position sorted.
         --threads=          [1] # threads to use for sorting.
         --sort_mem=         [1G] Amount of RAM to use for sorting.  
        --ref=              Reference fasta + index. 
        --output_dir=       Name the directory for output.  
         --output_prefix=    [filename] Prefix for the output file.
        --region=           chr1:100-200. Use to look @ reads only in this region. Highly recommended to use.
        --M_M=              <0|1> [0] 1= calculate Mapped_Mapped reads only. 
        --M_UM=             <0|1> [0] 1= calculate Mapped_UN-Mapped reads only. Works for LGT reads.
        --d=                [10000] Max Coverage per base. 
        --A=                <0|1> [1] 1= Count anomalous read pairs (LGT).
        --Qsub=             <0|1> [0] 1= Qsub the mpileup foreach input.\n";
}
#### Check input
if ( !$options{input} && $ARGV[0] ) { $options{input} = $ARGV[0]; }
if ( !$options{input} && !$options{input_list} ) { die "ERROR: Must give an input bam with --input=<BAM>. Try again.\n"; }
if ( $options{M_M}    && $options{M_UM} )        { die "ERROR: Can only use one: M_M or M_UM at once. Try again.\n"; }

## Setup options
my $qsub     = defined $options{Qsub}     ? $options{Qsub}       : "0";
my $sort     = defined $options{sort}     ? $options{sort}       : "0";
my $MM       = defined $options{MM}       ? $options{MM}         : "0";
my $MU       = defined $options{MU}       ? $options{MU}         : "0";
my $A        = defined $options{A}        ? $options{A}          : 1;
my $d        = defined $options{d}        ? $options{d}          : 100000;
my $region   = defined $options{region}   ? "'$options{region}'" : undef;
my $threads  = defined $options{threads}  ? "$options{threads}"  : 4;
my $sort_mem = defined $options{sort_mem} ? "$options{sort_mem}" : "1G";
my $ref      = defined $options{ref}      ? "$options{ref}"      : undef;
my $view    = "-hu";       ## Default
my $mpileup = "-Ad $d";    ## Default

if ( $MM == 1 ) {
    $view    = "-huF0xC";
    $mpileup = "-d $d";
}
if ( $MU == 1 )              { $view    = "-huf0x8 -F0x4"; }
if ( $A == 0 )               { $mpileup = "-d $d"; }
if ( defined $options{ref} ) { $mpileup = $mpileup . " -f $ref"; }

my $input = setup_input( \%options );

foreach my $bam (@$input) {
    my ( $filename, $directory ) = fileparse( $bam, ( '.srt.bam', '.bam' ) );
    my $prefix = $options{output_prefix} ? $options{output_prefix} : $filename;
    my $dir    = $options{output_dir}    ? $options{output_dir}    : $directory;
    if ( $dir =~ /\w+/ ) { run_cmd("mkdir -p -m u=rwx,g=rwx,o= $dir"); }
    my $out = ( $dir =~ /w+/ ) ? "$dir$prefix" : $prefix;
    my $cmd;
    if ( $sort == 1 && defined $options{region} ) { $cmd = "samtools sort -@ $threads -m $sort_mem -o $bam - | samtools view $view - $region | samtools mpileup $mpileup - > $out\.mpileup"; }
    elsif ( $sort == 1 && !$options{region} ) { $cmd = "samtools sort -@ $threads -m $sort_mem -o $bam - | samtools mpileup $mpileup - > $out\.mpileup"; }
    elsif ( defined $options{region} || $options{MM} == 1 || $options{MU} == 1 ) { $cmd = "samtools view $view $bam $region | samtools mpileup $mpileup - > $out\.mpileup"; }
    else                                                                         { $cmd = "samtools mpileup $mpileup $bam > $out\.mpileup"; }

    if ( $qsub == 1 ) { Qsub( { cmd => $cmd, threads => $threads, sub_mem => $sort_mem, wd => $dir, sub_name => "mpileup" } ); next; }
    else              { run_cmd($cmd); }
}
