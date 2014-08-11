#!/usr/local/bin/perl
use warnings;
use strict;
use LGTSeek;
use mk_dir;
use File::Basename;
use run_cmd;
use POSIX;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input=s', 'threads=s', 'output_prefix=s', 'output_dir=s', 'Qsub=s', 'sort_mem=s', 'sub_mem=s', 'sub_mail=s', 'help', )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) {
    die "Help: This script will take a bam to sort and index it.
      --input=                <BAM> to sort & index
      --output_prefix=        prefix.bam
      --output_dir=           directory to put output
      --sort_mem=        [1G] RAM / thread for sort. 
      --Qsub=            <0|1> [0] 1= Qsub the sort. 
        --sub_mem=       [5G] Needs to reflect changes in --sort_mem (sort_mem * threads + overhead).
        --threads=       [1] # of threads to use for sorting. 
        --sub_mail=      [0] 1= email user\@som.umaryland.edu when job is complete + with stats. Can also specify --sub_mail=specific\@email.foo
      --help\n";
}

if ( !$options{input} && !$ARGV[0] ) { die "Must pass an input bam to sort, use --input=<BAM> or pass \$ARGV[0]\n"; }
my $input = $options{input} ? $options{input} : $ARGV[0];
my $lgtseek = LGTSeek->new2( \%options );
my ( $bam, $path, $suf ) = fileparse( $input, @{ $lgtseek->{bam_suffix_list} } );
my $prefix = $options{output_prefix} ? $options{output_prefix} : $bam;
my $dir    = $options{output_dir}    ? $options{output_dir}    : $path;
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
        {   cmd      => "$sub",
            mem      => "$sort_mem",
            wd       => "$dir",
            sub_name => "sortBAM",
            threads  => "$threads",
            sub_mail => $options{sub_mail},
            no_gal   => 1,
        }
    );
    die "Job submitted to the grid.\n";
}

print STDERR "+++ Sorting bam +++\n";
my $cmd1 = "samtools sort -m $sort_mem -@ $threads $input $out\.psort";
run_cmd($cmd1);
print STDERR "+++ Indexing bam +++\n";
my $cmd2 = "samtools index $out\.psort.bam $out\.psort.bai";
run_cmd($cmd2);

