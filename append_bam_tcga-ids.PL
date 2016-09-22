#!/usr/bin/perl
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use warnings;
use strict;
use read_bam;
use mk_dir;

if ( !@ARGV ) { &help; }

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'output|O=s', 'output_dir|o=s', 'output_prefix|p=s', 'help|?' ) or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} ) { &help; }
if ( !defined $options{input} and !$ARGV[0] ) { die "Please pass an input bam to append TCGA-ids using --input=< /path/to.bam > or ARGV[0]\n"; }

my ( $fn, $path, $suf ) = fileparse( $input, ".bam" );
my $out_pref = defined $options{output_prefix} ? $options{output_prefix} : $fn;
my $out_dir  = defined $options{output_dir}    ? $options{output_dir}    : $path;
my $out      = defined $options{output}        ? $options{output}        : "$path$fn.bam";
my ( $out_final_fn, $out_final_dir, $out_final_suf ) = fileparse( $out, ".bam" );
mk_dir($out_final_dir);

my $input = defined $options{input} ? $options{input} : $ARGV[0];
my ( $header, $in ) = open_bam($input);
my $out = write_bam( $out, $header );

while (<$in>) {
    chomp( my $line = $_ );
    if ( $line =~ /AI:Z:([A-Za-z0-9\-]+)/ ) {
        my $a_id = $1;
        if ( defined $hash{$a_id} ) {
            print $out "$line\tPI:Z:$hash{$a_id}\n";
        }
        else {
            my $query = `cgquery "analysis_id=$a_id" | grep "participant_id"`;
            chomp($query);
            my @f = split( /:/, $query );
            $f[1] =~ s/\s+//g;
            $hash{$a_id} = $f[1];
            print $out "$line\tPI:Z:$hash{$a_id}\n";
        }
    }
    elsif ( $line =~ /FN:Z:([A-Za-z0-9\/\-\_\.]+)/ ) {
        my $file = $1;
        my $part_id;
        if    ( $file =~ /(SR[SR]{1}\d{6})/ )                 { $part_id = $1; }
        elsif ( $file =~ /devine-sample.*(COAD\d{1}_\d{1})/ ) { $part_id = $1; }
        elsif ( $file =~ /devine-sample.*(NORM\d{1}_\d{1})/ ) { $part_id = $1; }
        elsif ( $file =~ /german/ )                           { $part_id = "german"; }
        else                                                  { die "$line\n"; }
        print $out "$line\tPI:Z:$part_id\n";
    }
    else { die "Foo2: $line\n"; }
}

close $in;
close $out;
