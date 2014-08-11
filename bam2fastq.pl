#!/usr/local/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
use warnings;
no warnings 'uninitialized';
use strict;
use run_cmd;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input=s', 'name_sort=s', 'threads=i', 'sort_mem', 'fasta=s', 'no_nums=s', 'interleaved=s', 'output_dir=s', 'subdirs=s', 'Qsub=s', 'sub_mem=s', 'help|h', );
if ( $options{help} ) {
    die "This script will take a bam and output the fastq. 
    --input=		input bam or use ARGV[0]. Must use FULL file path names for Qsub.
      --name_sort=	<0|1> [0] 1=Sort the input bam by names first. 
      --sort_mem=	[5G] Amount of mem / thread to use for name sorting.
      --threads= 	[1] # of threads to use for name sorting.
    --fasta= 		<0|1> [0] 0=Fastq. 1=fasta
    --interleaved= 	<0|1> [0] 0=Seperate files for output(_1/_2). 1=Interleaved pairs.
    --no_nums=  	<0|1> [0] 0=Add _1 and _2 to the fastq. 1=Don't add numbers.
    --output_dir= 	[input dir]
      --subdirs=	<0|1> [0] 1= Make subdirectory to work from. 
    --output_list=	<0|1> [0] 1= Make list of output fastqs. Comma seperated, 1 pair/line. 
    --Qsub=			<0|1> [0] 1= Qsub. Must use FULL file path names for Qsub.
      --sub_mem=	[6G] Must reflect changes in --sort_mem.\n";
}
## Open Input
if ( !$ARGV[0] && !$options{input} ) { die "Must give an input bam with ARGV[0] or --bam=<input>.\n"; }

my $sub_mem = $options{sub_mem} ? $options{sub_mem} : "6G";
$options{input} = defined $options{input} ? $options{input} : $ARGV[0];
my $in = $options{input};
my ( $fn, $path, $suf ) = fileparse( $in, ( '.srt.bam', '_\w+_sort.bam', '.sort.bam', '.bam' ) );
my $out_dir = $options{output_dir} ? $options{output_dir} : $path;
my $subdirs = $options{subdirs}    ? $options{subdirs}    : 0;
if ( $options{subdirs} == 1 ) { $out_dir = "$out_dir/$fn"; }
run_cmd("mkdir -p $out_dir");

if ( $options{Qsub} == 1 ) {
    my $cmd = "/home/ksieber/scripts/bam2fastq.pl ";
    foreach my $key ( keys %options ) {
        next if ( $key =~ /Qsub/ );
        if ( $options{$key} ) { $cmd = $cmd . " --$key=$options{$key}" }
    }
    $fn =~ /(\w{1,10})$/;
    my $job_name = $1;
    Qsub(
        {   cmd     => $cmd,
            wd      => $out_dir,
            name    => $job_name,
            mem     => $options{sub_mem},
            threads => $threads,
            project => "jdhotopp-lab",
        }
    );
    die "Job submitted to the grid. Exiting.\n";
}

## Set Defaults
if ( !$options{no_nums} )     { $options{no_nums}     = 0; }
if ( !$options{fasta} )       { $options{fasta}       = 0; }
if ( !$options{interleaved} ) { $options{interleaved} = 0; }
if ( $options{no_nums} == 1 ) { $options{interleaved} = 0; }    ## If _1/_2 are NOT added, MUST prints reads to different files
my $name_sort_input = $options{name_sort}       ? $options{name_sort} : 0;
my $sort_mem        = $options{sort_mem}        ? $options{sort_mem}  : "5G";
my $threads         = defined $options{threads} ? $options{threads}   : "1";
if ( $name_sort_input == 1 ) {
    open( IN, "samtools sort -@ -m $sort_mem -no $in - | samtools view - |" ) || die "Can't open: $in because: $!\n";
}
else {
    open( IN, "-|", "samtools view $in" ) || die "Unable to open input $in.\n";
}

## Setup output
my $out1;
my $out2;
my $output = "$out_dir/$fn";
if ( $options{interleaved} == 0 ) {
    if ( $options{fasta} == 0 ) {
        open( $out1, ">", "$output\_1.fq" ) or die "Can't open output file: $output\_1.fq because: $!\n";
        open( $out2, ">", "$output\_2.fq" ) or die "Can't open output file: $output\_2.fq because: $!\n";
    }
    elsif ( $options{fasta} == 1 ) {
        open( $out1, ">", "$output\_1.fa" ) or die "Can't open output file: $output\_1.fa because: $!\n";
        open( $out2, ">", "$output\_2.fa" ) or die "Can't open output file: $output\_2.fa because: $!\n";
    }
}
elsif ( $options{interleaved} == 1 ) {
    if ( $options{fasta} == 0 ) {
        open( $out1, ">", "$output\.fq" ) or die "Can't open output file: $output\.fq because: $!\n";
        $out2 = $out1;
    }
    elsif ( $options{fasta} == 1 ) {
        open( $out1, ">", "$output\.fa" ) or die "Can't open output file: $output\.fa because: $!\n";
        $out2 = $out1;
    }
}

## Process .bam conversion to .fastq
my %foo;
while (<IN>) {
    chomp;
    next if $_ =~ /^@/;
    my @f = split;
    $foo{ $f[0] }++;
    if ( $foo{ $f[0] } == 2 ) {
        if ( $options{no_nums} == 0 ) { $f[0] =~ s/$f[0]/$f[0]\_2/; }
        _print( $f[0], $f[9], $f[10], $out2 );
        $foo{ $f[0] } = undef;
        next;
    }
    if ( $foo{ $f[0] } == 1 ) {
        if ( $options{no_nums} == 0 ) { $f[0] =~ s/$f[0]/$f[0]\_1/; }
        _print( $f[0], $f[9], $f[10], $out1 );
    }
}

print STDERR "Completed bam2fastq.pl on: $options{input}\n";

sub _print {
    my ( $read, $seq, $qual, $fh ) = @_;
    if ( $options{fasta} == 0 ) { print $fh "\@$read\n$seq\n\+\n$qual\n"; }
    if ( $options{fasta} == 1 ) { print $fh "\>$read\n$seq\n"; }
}

