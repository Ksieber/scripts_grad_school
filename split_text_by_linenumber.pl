#!/usr/bin/perl

=head1 NAME

split_text_by_linenumber.pl

=head1 SYNOPSIS

Splits input text into chunks of a certain size.

=head1 DESCRIPTION

Splits input text(s) into chunks of a certain size.

=head1 AUTHOR - Karsten Sieber

e-mail: ksieber@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

use warnings;
no warnings 'uninitialized';
use strict;
use Scalar::Util qw(reftype);
use File::Basename;
use run_cmd;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
my $results = GetOptions(
    \%options,       'input=s',        'line_count=s',  'Qsub=s',        'output_dir=s', 'subdirs=s',
    'output_list=s', 'samtools_bin=s', 'ergatis_dir=s', 'output_list=s', 'bin_dir=s',    'help|?'
) or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) {
    die "Help: This script will split bams into smaller bams based on the number of lines per file. 
      --input=         <TEXT>
      --line_count=    <lines per output> [1000]
      --Qsub=          <0|1> [0] 1= Qsub the split command.
      --output_dir=    Directory for output. 
      --subdirs        <0|1> [0] 1= Make a directory within the output_dir to place the output. 
      --output_list=   <0|1> [0] 1= Name the output with the list of files created.
      --help\n";
}

if ( !$options{input} ) { die "Must give an --input=<Text> to split." }
my $line_count = $options{line_count} ? $options{line_count} : "1000";
my $Qsub       = $options{Qsub}       ? $options{Qsub}       : 0;

my ( $fn, $path, $suf ) = fileparse( $options{input}, ( '.list', '.txt' ) );
if ( $path =~ /\.\// ) { chomp( $path = `pwd` ); $path = $path . "/"; }

my $output_dir = $options{output_dir} ? $options{output_dir} : $path;
my $subdirs    = $options{subdirs}    ? $options{subdirs}    : 0;
`mkdir -p $options{output_dir}`;

if ( $subdirs == 1 ) {
    $options{output_dir} = $options{output_dir} . "/$fn";
    `mkdir -p $options{output_dir}`;
}
if ( $Qsub == 1 ) {
    my $cmd = "/home/ksieber/scripts/split_text_by_linenumber.pl";
    foreach my $key ( keys %options ) {
        next if ( $key =~ /Qsub/ );
        if ( $options{$key} ) { $cmd = $cmd . " --$key=$options{$key}" }
    }
    $fn =~ /(\w{1,10})$/;    ## Grab the last 1-10 character of the input name to use as the job_name
    my $job_name = $1;
    Qsub(
        {   cmd  => $cmd,
            wd   => "$options{output_dir}",
            name => "$job_name",
        }
    );
    next;
}
else {
    my $count  = 0;
    my $i      = 0;
    my $output = "$options{output_dir}/$fn\_$count$suf";
    my @output_list;
    open( my $out, ">", "$output" )         || die "Can't open output: $output because: $!\n";
    open( my $in,  "<", "$options{input}" ) || die "Can't open input: $options{input} because: $!\n";
    push( @output_list, $output );
    while (<$in>) {

        if ( $i >= $line_count ) {
            close $out or die "Unable to close output: $output\n";
            $count += 1;
            $output = "$options{output_dir}/$fn\_$count$suf";
            open( $out, ">", "$output" ) || die "Can't open: $output because: $!\n";
            push( @output_list, $output );
            $i = 0;
        }
        print $out "$_";
        $i += 1;
    }
    close $out;
    close $in;

    if ( $options{output_list} == 1 ) {
        open( my $foo, ">", "$options{output_dir}/output.list" );
        print $foo join( "\n", @output_list );
        close $foo;
    }
}

__END__
