#!/usr/bin/perl -I /home/ksieber/scripts/ -I /home/ksieber/perl5/lib/perl5/
use warnings;
no warnings 'uninitialized';
use strict;
use File::Basename;
use run_cmd;
use mk_dir;
use print_call;
use setup_input;
use File::Basename;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'input_list|I=s', 'faidx=i', 'output_dir|o=s', 'output_prefix|p=s', 'output_list=s', 'Qsub|q=i', 'help|?', 'sub_name=s' )
    or die "Unrecognized comand line option. Please try again.\n";

if ( $options{help} ) {
    die "Help: This script will BWA index a fasta file for BWA aligner.
    --input|i=              Single Fasta file to index.
    --input_list=           Either a list of fasta files or a string of files comma seperated.
    --faidx=                <1|0> [1] 1= samtools faidx index fasta also. Needed for mpileup etc.
    --output_dir|o=         /path/for/output/ [/input/dir/]
    --output_prefix|p=              \$prefix for output names [input_prefix.fasta]
    --output_list=          /path/to/file.list to append ref output names to. 
    --Qsub|q=               1= Submit job to grid.
    --project=              [jdhotopp-lab] 
    --help|?
    \n";
}

if ( $options{Qsub} == 1 ) {
    $options{sub_name} = "Index";
    Qsub_script( \%options );
}

$options{faidx} = $options{faidx} ? $options{faidx} : 1;

my $inputs = setup_input( \%options );
foreach my $input (@$inputs) {
    my ( $fn, $dir, $suf ) = fileparse( $input, ( '.fa', '.fasta' ) );
    my $output_dir    = $options{output_dir}    ? $options{output_dir}    : $dir;
    my $output_prefix = $options{output_prefix} ? $options{output_prefix} : $fn;

    my $cmd = "bwa index -p $output_dir/$output_prefix$suf $input";
    run_cmd($cmd);

    if ( $options{faidx} == 1 ) {
        my $faidx = "samtools faidx $input";
        run_cmd($faidx);
    }

    if ( $options{output_list} ) {
        open( OUT, ">>", "$options{output_list}" ) or die "Error: Unable to open --output_list: $options{output_list}\n";
        print OUT "$output_dir/$output_prefix";
        close OUT;
    }
}

print_complete( \%options );
