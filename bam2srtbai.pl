#!/usr/local/bin/perl
use lib (
    '/home/ksieber/perl5/lib/perl5/',               '/home/ksieber/scripts/',
    '/local/projects-t3/HLGT/scripts/lgtseek/lib/', '/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/'
);
use warnings;
use strict;
use LGTSeek;
use mk_dir;
use File::Basename;
use run_cmd;
use POSIX;
use Cwd;

if ( !@ARGV ) { &help; }

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results
    = GetOptions( \%options, 'input|i=s', 'threads|t=s', 'output_prefix|p=s', 'output_dir|o=s', 'Qsub|q=s', 'sort_mem=s', 'sub_mem=s', 'sub_mail=s',
    'help|?', )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) { &help; }

if ( !$options{input} && !$ARGV[0] ) { die "Must pass an input bam to sort, use --input=<BAM> or pass \$ARGV[0]\n"; }
my $input = $options{input} ? $options{input} : $ARGV[0];
$options{input} = $input;
my $lgtseek = LGTSeek->new2( \%options );
my ( $bam, $path, $suf ) = fileparse( $input, @{ $lgtseek->{bam_suffix_list} } );
my $prefix = defined $options{output_prefix} ? $options{output_prefix} : "$bam\.psort";
if ( $path =~ /\.\// ) {
    my $cwd = getcwd;
    if ( -e "$cwd/$bam$suf" ) {
        $path  = $cwd;
        $input = "$cwd/$bam$suf";
    }
}
my $dir = $options{output_dir} ? $options{output_dir} : $path;
mk_dir($dir);
$dir =~ s/\/$//;
my $out      = "$dir/$prefix";
my $Qsub     = defined $options{Qsub} ? $options{Qsub} : 0;
my $threads  = defined $options{threads} ? $options{threads} : "1";
my $sort_mem = defined $options{sort_mem} ? $options{sort_mem} : "1G";
$sort_mem =~ /(\d+)G$/;
my $sub_mem = ( ceil( ( $1 * $threads ) * 1.1 ) + 1 ) . "G";

my $sub = "perl /home/ksieber/scripts/bam2srtbai.pl";
if ( $Qsub == 1 ) {
    foreach my $key ( keys %options ) {
        next if ( $key =~ /Qsub/ && !$options{input_list} );    ## If we are in the orignal call with input_list, we probably want to qsub each input
        if ( $options{$key} ) { $sub = $sub . " --$key=$options{$key}" }
        ;                                                       ## Build the command for all other options passed in @ original call
    }
    Qsub(
        {   cmd      => $sub,
            sub_mem  => $sub_mem,
            wd       => $dir,
            sub_name => "sortBAM",
            threads  => $threads,
            sub_mail => $options{sub_mail},
        }
    );
    die "Job submitted to the grid.\n";
}

## Sort the input into a unique tmp bam
my $random_int = ceil( rand(1) * 1000000 );
print STDERR "+++ Sorting bam +++\n";
run_cmd("samtools sort -m $sort_mem -@ $threads $input $out\_tmp_$random_int");
## Adjust header
print STDERR "+++ Adjusting header +++\n";
my $header = `samtools view -H $input`;
if   ( $header =~ /\@HD/ ) { $header =~ s/SO:\w+/SO:coordinate/; }
else                       { $header = "\@HD\tVN:1.1\tSO:coordinate\n" . $header; }
open( my $OUT, ">", "$out\_tmp_$random_int.header" ) or die "Error: Unable to open the output header file: $out\_tmp_$random_int.header\n";
print $OUT $header;
close $OUT;
run_cmd("samtools reheader $out\_tmp_$random_int\.header $out\_tmp_$random_int\.bam > $out\.bam");
run_cmd("/home/ksieber/bin/WC $out\.bam");
run_cmd("rm $out\_tmp_$random_int\.header $out\_tmp_$random_int\.bam");
## Index the bam
print STDERR "+++ Indexing bam +++\n";
run_cmd("samtools index $out\.bam");

sub help {
    die "\nHelp: This script will take a bam to sort and index it.
      --input|i=              <BAM> to sort & index
      --output_prefix|p=      <prefix>.bam
      --output_dir|o=         directory to put output
      --sort_mem=             [1G] RAM / thread for sort. 
      --threads|t=            [1] # of threads to use for sorting.
      --Qsub|q=               <0|1> [0] 1= Qsub the sort. 
        --sub_mem=            [5G] Needs to reflect changes in --sort_mem (sort_mem * threads + overhead). 
        --sub_mail=           [0] 1= email user\@som.umaryland.edu when job is complete + with stats. Can also specify --sub_mail=specific\@email.foo
      --help|?\n\n";
}
