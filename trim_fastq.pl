#!/usr/bin/perl
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use strict;
use warnings;
use mk_dir;
use File::Basename;
use Bio::SeqIO;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options;
my $results = GetOptions( \%options, 'input|i=s', 'range|r=s', 'help|?', 'output_dir|o=s', 'output_prefix|p=s', 'output|O=s' )
    or die "\n*** Error *** Unrecognized command line option. Please try again.\n\n";

if ( $options{help} )   { &help; }                                                                                       ## At the end of the script
if ( !$options{input} ) { die "Error: Must pass an input fasta file with --input=/some/file.fa Please try again.\n"; }
if ( !$options{range} ) { die "Error: Must pass an range to trim FOR in each FASTQ entry. ex: -r=1-50\n"; }

my ( $lower_range, $upper_range ) = ( split /\-/, $options{range} );
my $input = $options{input};
my ( $fn, $path, $suf ) = fileparse( $input, qr/\.[^\.]+/ );

my $output_dir    = defined $options{output_dir}    ? $options{output_dir}    : $path;
my $output_prefix = defined $options{output_prefix} ? $options{output_prefix} : $fn;
my $output        = defined $options{output}        ? $options{output}        : "$output_dir/$output_prefix\.fa";
my ( $o_fn, $o_path, $o_suf ) = fileparse( $output, qr/\.[^\.]+/ );
mk_dir($o_path);

my $OUT;
if ( defined $options{output_dir} or defined $options{output_prefix} or defined $options{output} ) {
    if ( -e $output ) { die "Error: Unsure what to name the output, $output already exists.\n"; }
    open( $OUT, ">$output" ) or die "Error: Unable to open the output.fa: $output\n";
}
else {
    $OUT = *STDOUT;
}

my $in         = Bio::SeqIO->new( -file => $input, -format => 'FASTQ' );
my $OUT_seq_fh = Bio::SeqIO->new( -fh   => \$OUT,  -format => 'FASTQ' );

while ( my $seq = $in->next_seq ) {
    my $new_seq = Bio::Seq::Quality->new(
        -id   => $seq->id,
        -seq  => $seq->subseq( $lower_range, $upper_range ),
        -qual => $seq->subqual_text( $lower_range, $upper_range ),
    );
    $OUT_seq_fh->write_seq($new_seq);
}
