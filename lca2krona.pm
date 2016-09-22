package lca2krona;
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use strict;
use warnings;
use mk_dir;
use File::Basename;
use Carp;
$Carp::MaxArgLen = 0;    ## Report full length error
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw( lca2krona );

=head2 &lca2krona
    Title       : lca2krona
    Usage       : my $krona = lca2krona({ lca => $lca_file, output_dir=> $dir, output_prefix => $prefix, remove_empty => 1 })
    Function    : Create an krona.html file for the input lca file.
    Returns     : Path to the krona.html file.
    Args        : \%hash
                    input
                    output_dir
                    output_prefix
                    remove_empty
                    col
=cut

sub lca2krona {
    my $options = shift;
    if ( !$options->{input} )        { confess "Error: &lca2krona didn't not receive an input properly.\n"; }
    if ( !( -e $options->{input} ) ) { confess "Error: &lca2krona input does not exist: $options->{input}\n"; }
    my ( $fn, $path, $suf ) = fileparse( $options->{input}, qr/\.[^\.]+/ );
    my $output_dir = defined $options->{output_dir} ? $options->{output_dir} : $path;
    mk_dir($output_dir);
    my $output_prefix = defined $options->{output_prefix} ? $options->{output_prefix} : $fn;
    my $col           = defined $options->{col}           ? $options->{col}           : "1";
    open( my $IN, "<", "$options->{input}" ) or confess "Error: &lca2krona unable to open input: $options->{input}\n";
    open( my $OUT, "| ktImportText -q -o $output_dir/$output_prefix\_krona.html - " ) or confess "Error: &lca2krona unable to open output: $output_dir/$output_prefix\_krona.html\n";

    while (<$IN>) {
        chomp;
        next if ( $_ =~ /^read_id/ );
        my @line = split( /\t/, $_ );
        my $original_lca = $line[$col];
        next if ( !$original_lca );
        my $new_lca = join( "\t", split( /;/, $original_lca ) );
        print $OUT "$new_lca\n";
    }
    close $IN;
    close $OUT;
    return "$output_dir/$output_prefix\_krona.html";
}

1;
