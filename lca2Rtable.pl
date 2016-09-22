#!/usr/bin/perl
use strict;
use warnings;
use lib ( '/home/ksieber/perl5/lib/perl5/', '/home/ksieber/scripts/' );
use read_in_list;
use File::Basename;
if ( !@ARGV ) { &help; }
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results
    = GetOptions( \%options, 'input|i=s', 'linker|L=s', 'lca_column|C=i', 'lca_cut|c=i', 'output|O=s', 'output_dir|0=s', 'output_prefix|p=s',
    'help|h|?' )
    or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} )            { &help; }
if ( !defined $options{input} )  { die "Error: Must use --input=<LCA file>\n"; }
if ( !defined $options{linker} ) { die "Error: Must use --linker=<file with Read-ids and participant/patient/analysis-ids\n"; }

my ( $in_fn, $in_dir, $in_suf ) = fileparse( $options{input}, qr/\.[^\.]+/ );
my $out_dir  = defined $options{output_dir}    ? $options{output_dir}    : $in_dir;
my $out_pref = defined $options{output_prefix} ? $options{output_prefix} : $in_fn;
my $output   = defined $options{output}        ? $options{output}        : "$out_dir$out_pref\_Rtable.txt";
my $OUT;
if ( defined $options{output} or defined $options{output_dir} or defined $options{output_prefix} ) {
    open( $OUT, ">", $output ) or die "Error: unable to open output: $output\n";
}
else { $OUT = *STDOUT; }

my $column  = defined $options{column}  ? $options{column}  : 1;
my $lca_cut = defined $options{lca_cut} ? $options{lca_cut} : 7;
my $links   = hash_in_data( $options{linker} );

my $lca_data;
my $otu_seen;
open( my $LCA, "<", $options{input} ) or die "Error: Couldn't open input: $options{input}\n";
while (<$LCA>) {
    chomp( my $line = $_ );
    my @split_line = split( /\t/, $_ );
    my $id  = $split_line[0];
    my $lca = $split_line[$column];
    next if ( !defined $lca or $lca !~ /\w+/ );
    my @split_lca = split( /;/, $lca );
    my @short_lca = ( $lca_cut eq "-1" ) ? @split_lca : splice( @split_lca, 0, $lca_cut );
    my $otu = $short_lca[-1];
    $otu =~ s/\s+/__/g;
    $otu =~ s/\//_/g;
    $otu =~ s/\-/_/g;
    $otu =~ s/\(//g;
    $otu =~ s/\)//g;
    $otu =~ s/\'//g;
    $otu =~ s/\"//g;
    $otu =~ s/\.//g;

    $lca_data->{ $links->{$id} }->{$otu}++;
    $otu_seen->{$otu}++;
}
close $LCA;

print $OUT "Ids";
foreach my $otu ( keys %{$otu_seen} ) {
    print $OUT "\t$otu";
}
print $OUT "\n";

foreach my $id ( keys %{$lca_data} ) {
    print $OUT "$id";
    foreach my $otu ( keys %{$otu_seen} ) {
        my $value = ( defined $lca_data->{$id}->{$otu} ) ? $lca_data->{$id}->{$otu} : "NA";
        print $OUT "\t$value";
    }
    print $OUT "\n";
}

sub help {
    die "This script will take an LGTSeek lca file and transform it into an Rtable.
	--input
	--Linker\n";
}
