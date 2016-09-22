#!/usr/bin/perl
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use strict;
use warnings;
use File::Basename;
use read_bam;
use Cwd;
use mk_dir;
use setup_input;
if ( !@ARGV ) { &help_full; }
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'output_dir|o=s', 'output_prefix|p=s', 'output|O=s', 'Qsub|q=i', 'help_full|?', 'min=i', 'max=i' )
    or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help_full} ) { &help_full; }                                             ## &help is @ the end of the script
if ( !$options{input} )    { die "Error: Must give an input bam with --input=.\n"; }
if ( !$options{min} )      { die "Error: Must give an input bam with --min=.\n"; }
if ( !$options{max} )      { die "Error: Must give an input bam with --max=\n"; }

my $input = $options{input};
my $min   = $options{min};
my $max   = $options{max};
my ( $in_fn, $in_dir, $in_suf ) = fileparse( $input, qr/\.[^\.]+/ );
my $out_dir  = ( defined $options{output_dir} )    ? $options{output_dir}    : $in_dir;
my $out_pref = ( defined $options{output_prefix} ) ? $options{output_prefix} : "$in_fn\_discord";
my $out_file = ( defined $options{output} )        ? $options{output}        : "$out_dir/$out_pref\_discord.bam";
my ( $final_out_fn, $final_out_dir, $final_out_suf ) = fileparse( $out_file, qr/\.[^\.]+/ );
mk_dir($final_out_dir);

my ( $header, $in ) = open_bam($input);
my $out = write_bam( $out_file, $header );
while ( my $read = read_bam($in) ) {
    my $print = 0;
    my $flag  = $read->{flag};

    # Parse reads MU or UM
    if ( !$flag->{qunmapped} and $flag->{munmapped} )  { $print = 1; goto PRINT; }
    if ( $flag->{qunmapped}  and !$flag->{munmapped} ) { $print = 1; goto PRINT; }

    # Parse reads RR, FF, RF
    if ( $flag->{qrev}  and $flag->{mrev} )  { $print = 1; goto PRINT; }
    if ( !$flag->{qrev} and !$flag->{mrev} ) { $print = 1; goto PRINT; }
    if ( $flag->{qrev}  and !$flag->{mrev} ) { $print = 1; goto PRINT; }

    # Parse reads on different chromosomes.
    if ( ( $read->{mate_ref} ne "=" ) and ( $read->{chr} ne $read->{mate_ref} ) ) { $print = 1; goto PRINT; }

    # Parse reads with unlikely distances
    if ( ( $read->{mate_ref} eq "=" ) or ( $read->{chr} eq $read->{mate_ref} ) ) {
        my $distance = abs( $read->{insert} );
        if ( $distance <= $min or $distance >= $max ) { $print = 1; goto PRINT; }
    }
PRINT: if ( $print == 1 ) {
        print $out "$read->{line}\n";
    }
}

close $in;
close $out;

sub help_full {
    die "This script will parse for discord reads.";
}
