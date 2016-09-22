#!/usr/bin/perl
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use strict;
use warnings;
use mk_dir;
use File::Basename;
use Bio::DB::Fasta;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options;
my $results = GetOptions( \%options, 'input|i=s', 'chr|c=s', 'range=s', 'help|?', 'output_dir|o=s', 'output_prefix|p=s', 'output|O=s' )
    or die "\n*** Error *** Unrecognized command line option. Please try again.\n\n";

if ( $options{help} )   { &help; }                                                                                        ## At the end of the script
if ( !$options{input} ) { die "Error: Must pass an input fasta file with --input=/some/file.fa Please try again.\n"; }
if ( !$options{chr} )   { die "Error: Must pass a \"chr\" to pull from the input.fa with --chr= Please try again.\n"; }

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

my $fasta_db = Bio::DB::Fasta->new($input);
my ( $lower_range, $upper_range ) = ( split /\-/, $options{range} ) if ( defined $options{range} );

if ( defined $options{chr} and $options{chr} =~ /all/ ) {
    if ( !defined $options{range} ) { die "Error: when --chr=\'all\' must use --range=\'#-#\'. Please try again. \n"; }
    my $stream = $fasta_db->get_PrimarySeq_stream;
    while ( my $seq = $stream->next_seq ) {
        print $OUT "\>" . $seq->id() . "\n" . $seq->subseq( $lower_range, $upper_range ) . "\n";
    }
}
else {
    my $new_seq = ( defined $options{range} ) ? $fasta_db->seq( $options{chr}, $lower_range => $upper_range ) : $fasta_db->seq( $options{chr} );
    print $OUT "\>$options{chr}";
    print $OUT ":$lower_range\-$upper_range" if ( defined $options{range} );
    print $OUT "\n";
    print $OUT "$new_seq\n";
}

close $OUT;

sub help {
    die "\n	This script will pull a fasta entry and/or range from a fasta file and 
	print it to a new fasta if any --output* is given else it prints to STDOUT.

	--input|i=				/some/input.fa to pull the sequence from.
	--chr|c=				Desired Fasta entry \"chr\" from the input file.
	  ** If --chr='all' script will pull --range for all chr.
	--range=				<###-###> Optional positional range for the chr. 
	--output|O=				Exact output desired. ex: /foo/bar/out.fa
	--output_dir|o=				Directory for output.
	--output_prefix|p=			Prefix for output.
	--help\n\n";
}
