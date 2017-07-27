#!/usr/bin/perl

=head1 NAME

split_bam_by_linenumber.pl

=head1 SYNOPSIS

Splits input bam(s) into chunks of a certain size.

=head1 DESCRIPTION

Splits input bam(s) into chunks of a certain size.

=head1 AUTHOR - Karsten Sieber

e-mail: ksieber@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut
use warnings;
use strict;
use local::lib;
use Scalar::Util qw(reftype);
use File::Basename;
use POSIX;
use run_cmd;
use mk_dir;
use linecount;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
my $results
    = GetOptions( \%options, 'input|i=s', 'lists|l=i', 'Qsub|q=s', 'output_dir|o=s', 'output_prefix|p=s', 'subdirs=s', 'sub_mem=s', 'samtools_bin=s', 'ergatis_dir=s', 'output_list=s', 'bin_dir=s',
    'help|?' )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) {
    die "Help: This script will split text randomly into the specified number of lists. 
      --input=          <Text List>
      --lists=          Number of lists to create.
      --Qsub=           <0|1> [0] 1= Qsub the split command.
      --output_dir=     Directory for output.
      --output_prefix= 
      --subdirs         <0|1> [0] 1= Make a directory within the output_dir to place the output. 
      --output_list=    <0|1> [0] 1= Name the output with the list of files created.
      --help\n";
}

if ( !$options{input} ) { die "Must give an --input=<Text> to split." }
if ( !$options{lists} ) { die "Must give an --lists=<#> to split into." }
$options{sub_mem} = $options{sub_mem} ? $options{sub_mem} : "100M";
if ( $options{Qsub} ) { Qsub_script( \%options ); }

my ( $fn, $path, $suf ) = fileparse( $options{input}, ( '.list', '.txt' ) );
if ( $path =~ /\.\// ) { chomp( $path = `pwd` ); $path = $path . "/"; }

my $output_dir    = defined $options{output_dir}    ? $options{output_dir}    : $path;
my $output_prefix = defined $options{output_prefix} ? $options{output_prefix} : $fn;
my $subdirs       = defined $options{subdirs}       ? $options{subdirs}       : 0;
if ( $subdirs == 1 ) { $output_dir = "$output_dir/$output_prefix"; }
mk_dir("$output_dir");

my $number_of_lists = $options{lists};

my $input_linecount = wc( $options{input} );
my $max_prints      = ceil( $input_linecount / $number_of_lists );

open( my $IN, "<", "$options{input}" ) or die "Error: Unable to open input: $options{input} because: $!\\n";

my $fh = {};
for ( my $i = 0; $i < $number_of_lists; $i++ ) {
    open( $fh->{$i}->{fh}, ">", "$output_dir/$output_prefix\_$i.list" ) or die "Error: Unable to open output_list: $output_dir/$output_prefix\_$i.list because: $!\n";
    $fh->{$i}->{open} = 1;
}

while (<$IN>) {
    my $random_fh = int( rand($number_of_lists) );
    if ( $random_fh >= $number_of_lists ) { $random_fh = 0; }
    my $to_print = 1;
    while ( $to_print == 1 ) {
        if ( !$fh->{$random_fh}->{open} ) { print "fail: $random_fh\n"; }
        if ( $fh->{$random_fh}->{open} == 1 ) {
            print { $fh->{"$random_fh"}->{fh} } $_;
            $to_print = 0;
            next;
        }
        elsif ( $random_fh >= $number_of_lists ) {
            $random_fh = 0;
        }
        elsif ( $random_fh >= 0 ) {
            $random_fh++;
            if ( $random_fh >= $number_of_lists ) {
                $random_fh = 0;
            }
        }
    }

    # Update fh status
    $fh->{$random_fh}->{prints}++;
    if ( $fh->{$random_fh}->{prints} == $max_prints ) {
        $fh->{$random_fh}->{open} = 2;
    }
}

close $IN;
for ( my $i = 0; $i < $number_of_lists; $i++ ) {
    close $fh->{$i}->{fh};
}

