#!/usr/local/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
use warnings;
no warnings 'uninitialized';
use strict;
use run_cmd;
use mk_dir;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,        'input|i=s',         'name_sort=s', 'threads|t=i', 'sort_mem=s', 'fasta=i', 'no_nums=i', 'interleaved=i',
    'output_dir|o=s', 'output_prefix|p=s', 'subdirs=s',   'Qsub|q=s',    'sub_mem=s',  'help|?',
) or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} ) {
    die "
    This script will take a bam and output the fastq.
    _____________________________________________________________________________________________ 
    --input|i=                <input.bam>, or use ARGV[0], or read STDIN. Must use FULL file path names for Qsub.
      --name_sort=            <0|1> [0] 1=Sort the input bam by names first. Must be name sorted otherwise.
      --sort_mem=             [1G] Amount of mem / thread to use for name sorting.
      --threads|t=              [1] # of threads to use for name sorting.
    _____________________________________________________________________________________________
    --fasta=                <0|1> [0] 0=Fastq. 1=fasta
    --interleaved=          <0|1> [1] 0=Seperate files for output(_1/_2). 1=Interleaved pairs.
    --no_nums=              <0|1> [0] 0=Add _1 and _2 to the fastq. 1=Don't add numbers.
      **Highly suggested to add _1/_2 when output is interleaved. **
    _____________________________________________________________________________________________
    --output_dir|o=           </path/for/output/> [input dir]
      --output_prefix|p=        <prefix.fq> [*prefix*.input]
      --subdirs=              <0|1> [0] 1= Make subdirectory to work from. 
    _____________________________________________________________________________________________ 
    --Qsub|q=                 <0|1> [0] 1= Qsub. Must use FULL file path names for Qsub.
      --sub_mem=              [1G] Must reflect changes in --sort_mem.
    --help|?
    _____________________________________________________________________________________________\n";
}
## Open Input
if ( !$ARGV[0] && !$options{input} && $options{Qsub} == 1 )      { die "Must give an input bam with ARGV[0] or --bam=<input.bam> when submitting to the grid.\n"; }
if ( !$ARGV[0] && !$options{input} && $options{name_sort} == 1 ) { die "Must give an input bam with ARGV[0] or --bam=<input.bam> when name_sort=1.\n"; }

my $sub_mem = $options{sub_mem} ? $options{sub_mem} : "1G";
$options{input} = defined $options{input} ? $options{input} : $ARGV[0];
my $in = ( -e $options{input} ) ? $options{input} : "<STDIN>";
my ( $fn, $path, $suf ) = ( -e $in ) ? fileparse( $in, ( '.srt.bam', '_\w+_sort.bam', '.sort.bam', '.bam' ) ) : ( "bam2fastq", "./", ".fa" );
my $out_dir = defined $options{output_dir}    ? $options{output_dir}    : $path;
my $prefix  = defined $options{output_prefix} ? $options{output_prefix} : $fn;
my $subdirs = defined $options{subdirs}       ? $options{subdirs}       : 0;
my $threads = defined $options{threads}       ? $options{threads}       : "1";
my $name_sort_input = ( -e $in && $options{name_sort} ) ? $options{name_sort} : 0;
my $sort_mem = defined $options{sort_mem} ? $options{sort_mem} : "1G";

if ( $options{subdirs} == 1 ) { $out_dir = "$out_dir/$fn"; }
mk_dir($out_dir);

if ( $options{Qsub} == 1 ) {
    my $cmd = "/home/ksieber/scripts/bam2fastq.pl ";
    foreach my $key ( keys %options ) {
        next if ( $key =~ /Qsub/ );
        if ( $options{$key} ) { $cmd = $cmd . " --$key=$options{$key}" }
    }
    Qsub(
        {   cmd      => $cmd,
            wd       => $out_dir,
            sub_name => "bam2fastq",
            sub_mem  => $options{sub_mem},
            threads  => $threads,
            project  => "jdhotopp-lab",
        }
    );
    die "Job submitted to the grid. Exiting.\n";
}

## Set Defaults
if ( !$options{no_nums} ) { $options{no_nums} = 0; }
if ( !$options{fasta} )   { $options{fasta}   = 0; }
if ( $options{interleaved} != 0 && $options{interleaved} != 1 ) { print STDERR "FOO\n"; $options{interleaved} = 1; }

my $IN;
if ( -e $in ) {
    if ( $name_sort_input == 1 ) {
        open( $IN, "-|", "samtools sort -@ $threads -m $sort_mem -no $in - | samtools view -" ) or die "Can't open: $in because: $!\n";
    }
    else {
        open( $IN, "-|", "samtools view $in" ) || die "Unable to open input $in.\n";
    }
}
else {
    $IN = *STDIN;
}

## Setup output
my $out1;
my $out2;
my $output = "$out_dir/$prefix";
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
while (<$IN>) {
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
close $IN;

print STDERR "*** Completed bam2fastq.pl on: $in ***\n";

sub _print {
    my ( $read, $seq, $qual, $fh ) = @_;
    if ( $options{fasta} == 0 ) { print $fh "\@$read\n$seq\n\+\n$qual\n"; }
    if ( $options{fasta} == 1 ) { print $fh "\>$read\n$seq\n"; }
}

