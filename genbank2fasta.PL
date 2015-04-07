#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
use warnings;
use strict;
use lib ( '/home/ksieber/perl5/lib/perl5/', '/local/projects-t3/HLGT/scripts/lgtseek/lib/', '/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/' );

=head1 NAME

genbank2fasta.PL

=head1 SYNOPSIS

Create fasta file from a genbank file.

=head1 EXAMPLE

genbank2fasta.PL important.gb
genbank2fasta.PL -i important.gb --output=/full/path/name/output.fa

=head1 AUTHOR - Karsten Sieber

e-mail: Karsten.sieber@gmail.com

=cut

use Bio::SeqIO;
use File::Basename;
use Cwd;
use mk_dir;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'output|O=s', 'output_dir|s', 'output_prefix=s', 'help|?' ) or die "\n*** Error *** Unrecognized command line option. Please try again.\n\n";

if ( !$options{input} and !@ARGV ) { die "Please pass the script an input genbank file using --input|i= or ARGV.\n\n"; }
my $genbank = defined $options{input} ? $options{input} : $ARGV[0];
if ( !-e $genbank ) { die "***ERROR*** The input.genbank file doens't exist: $genbank\n\n"; }

my ( $fn, $path, $suffix ) = fileparse( $genbank, qr/\.[^\.]+/ );
my $out_prefix = defined $options{output_prefix} ? $options{output_prefix} : $fn;
my $out_dir    = defined $options{output_dir}    ? $options{output_dir}    : $path;
if ( $out_dir =~ /^\.\/$/ ) {
    my $cwd = getcwd;
    if ( -e "$cwd/$fn$suffix" and !-e "$cwd/$fn\.fa" ) { $out_dir = $cwd; }
}
mk_dir($out_dir);
my $out_name = "$out_dir/$fn\.fa";
my $out      = Bio::SeqIO->new( -file => ">$out_name", -format => 'fasta' );
my $in       = Bio::SeqIO->new( -file => $genbank, -format => 'genbank' );

while ( my $seq = $in->next_seq() ) { $out->write_seq($seq); }
