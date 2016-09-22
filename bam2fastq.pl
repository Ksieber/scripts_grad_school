#!/usr/local/bin/perl 
use lib ( 
  '/home/ksieber/perl5/lib/perl5/', 
  '/home/ksieber/scripts/', 
  '/local/projects-t3/HLGT/scripts/lgtseek/lib/',
  '/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/'
  );
use warnings;
no warnings 'uninitialized';
use strict;
use Cwd;
use run_cmd;
use read_in_list;
use mk_dir;
use Bio::Util::DNA qw(:all);    ## 05.20.15
use parse_flag;                 ## 05.20.15
use File::Basename;

if ( -t STDIN and not @ARGV ) { &help; }

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,    'input|i=s', 'name_sort=s', 'threads|t=i', 'sort_mem=s', 'fasta=i',     'no_nums=i',  'interleaved=i', 'output_dir|o=s', 'output_prefix|p=s',
    'output|O=s', 'subdirs=s', 'Qsub|q=s',    'sub_mem=s',   'help|?',     'good_list=s', 'bad_list=s', 'gzip=i'
) or die "Unrecognized command line option. Please try agian.\n";

## Open Input
if ( $options{help} ) { &help; }
if ( !$ARGV[0] && !$options{input} && $options{Qsub} == 1 ) {
    die "Must give an input bam with ARGV[0] or --bam=<input.bam> when submitting to the grid.\n";
}
if ( !$ARGV[0] && !$options{input} && $options{name_sort} == 1 ) {
    die "Must give an input bam with ARGV[0] or --bam=<input.bam> when name_sort=1.\n";
}

## Read in lists of good or bad read-ids
my $good_ids = &hash_in_list( $options{good_list} ) if ( defined $options{good_list} );
my $bad_ids  = &hash_in_list( $options{bad_list} )  if ( defined $options{bad_list} );

my $sub_mem = $options{sub_mem} ? $options{sub_mem} : "1G";
$options{input} = defined $options{input} ? $options{input} : $ARGV[0];
my $in = ( -e $options{input} ) ? $options{input} : "<STDIN>";
my ( $fn, $path, $suf ) = ( -e $in ) ? fileparse( $in, ( '.psort.bam', '.srt.bam', '_\w+_sort.bam', '.sort.bam', '.bam' ) ) : ( "bam2fastq", "./", ".fa" );
if ( $path !~ /\w+/ ) { $path = getcwd; }
my $out_dir = defined $options{output_dir}    ? $options{output_dir}    : $path;
my $prefix  = defined $options{output_prefix} ? $options{output_prefix} : $fn;
my $output  = defined $options{output}        ? $options{output}        : "$out_dir/$prefix";
my $subdirs = defined $options{subdirs}       ? $options{subdirs}       : 0;
my $threads = defined $options{threads}       ? $options{threads}       : 1;
my $gzip    = defined $options{gzip}          ? $options{gzip}          : 0;
my $name_sort_input = ( -e $in && $options{name_sort} ) ? $options{name_sort} : 0;
my $sort_mem = defined $options{sort_mem} ? $options{sort_mem} : "1G";

if ( $options{subdirs} == 1 ) { $out_dir = "$out_dir/$fn"; }
mk_dir($out_dir);

if ( $options{Qsub} == 1 ) {
    my $cmd = "/home/ksieber/scripts/bam2fastq.pl ";
    foreach my $key ( keys %options ) {
        next if ( $key =~ /Qsub/ );
        if ( defined $options{$key} ) { $cmd = $cmd . " --$key=$options{$key}" }
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
if ( not $options{output} and not $options{output_dir} and not $options{output_prefix} and not $options{interleaved} ) { $options{interleaved} = 1; }
if ( !defined $options{no_nums} )     { $options{no_nums}     = 1; }
if ( !defined $options{fasta} )       { $options{fasta}       = 0; }
if ( !defined $options{interleaved} ) { $options{interleaved} = 0; }

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
if ( ( not $options{output} and not $options{output_dir} and not $options{output_prefix} ) and ( $options{interleaved} == 1 ) ) {
    $out1 = *STDOUT;
    $out2 = *STDOUT;
}
elsif ( $options{interleaved} == 0 ) {
    if ( $options{fasta} == 0 ) {
        if ( defined $gzip and $gzip == 1 ) {
            open( $out1, "| gzip -c - > $output\_1.fq.gz" ) or die "Can't open output file: $output\_1.fq.gz because: $!\n";
            open( $out2, "| gzip -c - > $output\_2.fq.gz" ) or die "Can't open output file: $output\_2.fq.gz because: $!\n";
        }
        else {
            open( $out1, ">", "$output\_1.fq" ) or die "Can't open output file: $output\_1.fq because: $!\n";
            open( $out2, ">", "$output\_2.fq" ) or die "Can't open output file: $output\_2.fq because: $!\n";
        }
    }
    elsif ( $options{fasta} == 1 ) {
        if ( defined $gzip and $gzip == 1 ) {
            open( $out1, "| gzip -c - > $output\_1.fa.gz" ) or die "Can't open output file: $output\_1.fa.gz because: $!\n";
            open( $out2, "| gzip -c - > $output\_2.fa.gz" ) or die "Can't open output file: $output\_2.fa.gz because: $!\n";
        }
        else {
            open( $out1, ">", "$output\_1.fa" ) or die "Can't open output file: $output\_1.fa because: $!\n";
            open( $out2, ">", "$output\_2.fa" ) or die "Can't open output file: $output\_2.fa because: $!\n";
        }
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
    my @f        = split;
    my $query_id = $f[0];
    if ( defined $options{good_list} and !$good_ids->{$query_id} ) { next; }
    if ( defined $options{bad_list}  and $bad_ids->{$query_id} )   { next; }

    # Reverse complement the sequence and quality score if needed
    my $sequence;
    my $quality_score;
    my $flag = parse_flag( $f[1] );
    if ( $flag->{qrev} ) {
        $sequence      = ${ reverse_complement( \$f[9] ) };
        $quality_score = reverse( $f[10] );
    }
    else {
        $sequence      = $f[9];
        $quality_score = $f[10];
    }

    # Print the output
    if ( $flag->{paired} ) {
        if ( $flag->{first} ) {
            if ( defined $options{no_nums} and $options{no_nums} == 0 ) { $f[0] =~ s/$f[0]/$f[0]\_2/; }
            _print( $f[0], $sequence, $quality_score, $out1 );
        }
        elsif ( !$flag->{first} ) {
            if ( defined $options{no_nums} and $options{no_nums} == 0 ) { $f[0] =~ s/$f[0]/$f[0]\_1/; }
            _print( $f[0], $sequence, $quality_score, $out2 );
        }
    }
    else {
        if ( defined $options{no_nums} and $options{no_nums} == 0 ) { $f[0] =~ s/$f[0]/$f[0]\_2/; }
        _print( $f[0], $sequence, $quality_score, $out1 );
    }

}
close $IN;

print STDERR "*** Completed bam2fastq.pl on: $in ***\n";

sub _print {
    my ( $read, $seq, $qual, $fh ) = @_;
    if ( $options{fasta} == 0 ) { print $fh "\@$read\n$seq\n\+\n$qual\n"; }
    if ( $options{fasta} == 1 ) { print $fh "\>$read\n$seq\n"; }
}

sub help {
    die "
    This script will take a bam and output the fastq.
    _____________________________________________________________________________________________ 
    --input|i=                  <input.bam>, or use ARGV[0], or read STDIN. Must use FULL file path names for Qsub.
      --name_sort=              <0|1> [0] 1=Sort the input bam by names first. Must be name sorted otherwise.
      --sort_mem=               [1G] Amount of mem / thread to use for name sorting.
      --threads|t=              [1] # of threads to use for name sorting.
    --good_list=                < /path/to/listof/desired/read-ids > 
    --bad_list=                 < /path/to/listof/UNdesired/read-ids > 
    _____________________________________________________________________________________________
    --fasta=                    <0|1> [0] 0=Fastq. 1=fasta
    --gzip=                     <0|1> [0] 1=gzip output
    --interleaved=              <0|1> [0] 0=Seperate files for output(_1/_2). 1=Interleaved pairs.
    --no_nums=                  <0|1> [1] 0=Add _1 and _2 to the fastq. 1=Don't add numbers.
      **Highly suggested to add _1/_2 when output is interleaved. **
    _____________________________________________________________________________________________
    If interleaved=1, output goes to STDOUT by default for piping
    --output|O=                 /full/path/name.out  
      --output_dir|o=           </path/for/output/> [input dir]
      --output_prefix|p=        <prefix.fq> [*prefix*.input]
      --subdirs=                <0|1> [0] 1= Make subdirectory to work from. 
    _____________________________________________________________________________________________ 
    --Qsub|q=                 <0|1> [0] 1= Qsub. Must use FULL file path names for Qsub.
      --sub_mem=              [1G] Must reflect changes in --sort_mem.
    --help|?
    _____________________________________________________________________________________________\n";
}
