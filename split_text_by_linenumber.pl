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
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use warnings;
no warnings 'uninitialized';
use strict;
use Scalar::Util qw(reftype);
use File::Basename;
use run_cmd;
use mk_dir;
use Cwd;

if ( !@ARGV ) { &help; }

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
my $results = GetOptions(
    \%options,       'input|i=s',      'line_count=s',  'Qsub|q=s',  'output_dir|o=s', 'output_prefix|p=s', 'output_suffix=s', 'subdirs=s',
    'output_list=i', 'samtools_bin=s', 'ergatis_dir=s', 'bin_dir=s', 'help|?',         'gzip|gz=i',
) or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) { &help; }
if ( !$options{input} and -t STDIN ) { die "Error: Must give an input via --input=</path/to/input.txt> or from STDIN.\n"; }
if ( !$options{input} and !$options{output_prefix} ) { die "Error: Must give an --output_prefix= when streaming STDIN.\n"; }

my $line_count = $options{line_count} ? $options{line_count} : "1000";
my $Qsub       = $options{Qsub}       ? $options{Qsub}       : "0";

my $tmp_out_fn;
my $tmp_out_dir;
my $tmp_out_suf;
if ( defined $options{input} ) { ( $tmp_out_fn, $tmp_out_dir, $tmp_out_suf ) = fileparse( $options{input}, qr/\.[^\.]+/ ); }
else                           { $tmp_out_fn = "split"; $tmp_out_dir = getcwd; $tmp_out_suf = ".txt"; }

my $gzip     = defined $options{gzip}          ? $options{gzip}          : "0";
my $out_dir  = defined $options{output_dir}    ? $options{output_dir}    : $tmp_out_dir;
my $out_pref = defined $options{output_prefix} ? $options{output_prefix} : $tmp_out_fn;
my $out_suf  = defined $options{output_suffix} ? $options{output_suffix} : $tmp_out_suf;
if ( $gzip == 1 and $out_suf !~ /\.gz$/ ) { $out_suf = $out_suf . ".gz"; }
my $subdirs = $options{subdirs} ? $options{subdirs} : "0";
if ( $subdirs == 1 ) { $out_dir = $out_dir . $out_pref; }
mk_dir($out_dir);

if ( $Qsub == 1 ) {
    $options{output_dir} = $out_dir;
    my $cmd = "/home/ksieber/scripts/split_text_by_linenumber.pl";
    foreach my $key ( keys %options ) {
        next if ( $key =~ /Qsub/ );
        if ( $options{$key} ) { $cmd = $cmd . " --$key=$options{$key}" }
    }
    $out_pref =~ /(\w{1,10})$/;    ## Grab the last 1-10 character of the input name to use as the job_name
    my $job_name = $1;
    my $threads = $gzip == 1 ? 2 : 1;
    Qsub(
        {   cmd      => $cmd,
            wd       => $out_dir,
            sub_name => "splitTXT",
            threads  => $threads,
        }
    );
    next;
}
else {
    my $count  = 0;
    my $i      = 0;
    my $output = "$out_dir/$out_pref\_$count$out_suf";
    my @output_list;
    my $out;
    if ( $gzip == 1 ) { open( $out, "| gzip -c - > $output" ) || die "Can't open output: $output because: $!\n"; }
    else              { open( $out, ">", "$output" ) || die "Can't open output: $output because: $!\n"; }

    my $in;
    if ( defined $options{input} and -e $options{input} ) {
        open( $in, "<", "$options{input}" ) || die "Can't open input: $options{input} because: $!\n";
    }
    elsif ( !-t STDIN ) {
        $in = *STDIN;
    }
    push( @output_list, $output );
    while (<$in>) {

        if ( $i >= $line_count ) {
            close $out or die "Unable to close output: $output\n";
            $count += 1;
            $output = "$out_dir/$out_pref\_$count$out_suf";
            if ( $gzip == 1 ) { open( $out, "| gzip -c - > $output" ) || die "Can't open output: $output because: $!\n"; }
            else              { open( $out, ">", "$output" ) || die "Can't open output: $output because: $!\n"; }
            push( @output_list, $output );
            $i = 0;
        }
        print $out "$_";
        $i += 1;
    }
    close $out;
    close $in;

    if ( $options{output_list} == 1 ) {
        open( my $foo, ">", "$out_dir/output.list" );
        print $foo join( "\n", @output_list );
        close $foo;
    }
}

sub help {
    die "\nHelp: This script will split text into smaller bams based on the number of lines per file. 
      --input|i=        <TEXT> or through STDIN
      --line_count=     <lines per output> [1000]
      --Qsub=           <0|1> [0] 1= Qsub the split command.
      --output_dir=     Directory for output. 
      --subdirs         <0|1> [0] 1= Make a directory within the output_dir to place the output. 
      --output_prefix=  \${prefix}\_#.\${suffix} 
      --output_suffix=  \${prefix}\_#.\${suffix}
      --gzip|gz=        <0|1> [0] 1= gzip each output. 
      --output_list=    <0|1> [0] 1= Name the output with the list of files created.
      --help\n\n";
}
