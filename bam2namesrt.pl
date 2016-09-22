#!/usr/bin/perl -I /home/ksieber/scripts/ -I /home/ksieber/perl5/lib/perl5/
use warnings;
use strict;
use File::Basename;
use setup_input;
use run_cmd;
use print_call;
use POSIX;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'input_list|I=s', 'output_prefix=s', 'output_dir|o=s', 'threads|t=i', 'sort_mem=s', 'sub_mem=s', 'Qsub|q=i', 'help|?', 'sub_mail=s', )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) {
    die "Help: This script will take a bam to sort based on names. 
      --input|i=           <BAM> to name sort.
      --input_list|I=      Either a file with a list 1 bam/line or comma seperated string of bams. 
      --output_dir|o=      directory to put output
      --output_prefix=     prefix.bam
      --sort_mem=          [1G] RAM / thread for sort.
      --threads|t=         [1] # of threads to use for sorting.  
      --Qsub|q=            <0|1> [0] 1= Qsub the sort. 
        --sub_mem=         [5G] Needs to reflect changes in --sort_mem (sort_mem * threads + overhead).
        --sub_mail=        [0] 1= email user\@som.umaryland.edu when job is complete + with stats. Can also specify --sub_mail=specific\@email.foo
      --help\n";
}

$options{input} = defined $options{input} ? $options{input} : $ARGV[0];
if ( !$options{input} && !$options{input_list} ) { die "Error: Must pass an input file with --input, --input_list, or $ARGV[0]. Please try again.\n"; }

if ( $options{Qsub} ) {
    $options{sub_mem}  = defined $options{sub_mem}  ? "$options{sub_mem}"  : "5G";
    $options{sub_name} = defined $options{sub_name} ? "$options{sub_name}" : "NameSortBam";
    my $original_sub_mem;
    my $original_sort_mem;
    if ( $options{sub_mem} ) {
        if ( $options{sub_mem} =~ /^(\d+)[K|M|G]$/ ) { $original_sub_mem = $1; }
    }
    if ( $options{sort_mem} ) {
        if ( $options{sort_mem} =~ /^(\d+)[K|M|G]$/ ) { $original_sort_mem = $1; }
    }
    if ( $original_sub_mem < ( $original_sort_mem * $options{threads} ) ) {
        $options{sub_mem} = ( ceil( ( $original_sort_mem * $options{threads} ) * 1.1 ) ) + 1 . "1G";
    }
    Qsub_script( \%options );
}

my $input_list = setup_input( \%options );
foreach my $input (@$input_list) {
    my ( $bam, $path, $suf ) = fileparse( $input, ( '-psrt\.bam', '\-psort\.bam', '\_psort\.bam', '\_pos-sort\.bam', '\.srt\.bam', '\.sort\.bam', '\.bam' ) );
    my $prefix = $options{output_prefix} ? $options{output_prefix} : $bam;
    my $dir    = $options{output_dir}    ? $options{output_dir}    : $path;
    $dir =~ s/\/$//;
    my $out      = "$dir/$prefix";
    my $threads  = defined $options{threads} ? "$options{threads}" : "1";
    my $sort_mem = $options{sort_mem} ? $options{sort_mem} : "1G";

    my $cmd = "samtools sort -n -@ $threads -m $sort_mem $input $out\_name-sort";
    run_cmd($cmd);
}

print_complete( \%options );
