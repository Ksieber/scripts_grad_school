#!/usr/bin/perl -I /home/ksieber/scripts/ -I /home/ksieber/perl5/lib/perl5/

=head1 NAME

get_genbank_entries.pl

=head1 SYNOPSIS

USAGE:  get_genbank_entries.pl --input=NC_0000001

=head2 OPTIONS

=over 24

=item B<--input|acc|i=>

Single accession

=item B<--input_list|list|I=>

List of accessions, either comma seperate or a file with 1 per line.

=item B<--type=>

"FASTA" or "gb" file format.

=item B<--output=>

/full/path/and/name/for/output.fa

=item B<--output_dir|o=>

/full/path/for/output/

=item B<--output_prefix|p=>

${output_prefix}.fa

=back

=head1 DESCRIPTION

get_genbank_entries.pl - Retrieves genbank records from genbank.

=cut

use strict;
use warnings;
use Bio::DB::EUtilities;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use setup_input;
use Cwd;

# Process the options list.
my %options;
my $results = GetOptions( \%options, 'input|acc|i=s', 'input_list|list|I=s', 'type=s', 'batch', 'help|h|?', 'output|O=s', 'output_dir|o=s', 'output_prefix|p=s', ) or &pod2usage();

&pod2usage( -exitval => 1, -verbose => 2, -output => \*STDOUT ) if ( $options{help} or !( !$options{input} or !$options{input_list} ) );

my $rettype = defined $options{type} ? $options{type} : 'FASTA';
$rettype =~ s/fasta/FASTA/;    ## force caps for EUtils
my $out_dir = defined $options{output_dir} ? $options{output_dir} : getcwd;
$out_dir =~ s/\/{2,}/\//g;
$out_dir =~ s/\/$//;
my $out_suffix = ( $rettype eq 'FASTA' ) ? ".fa" : ".gbk";

# Loop through the IDs and pull down the genbank file.
# Write the results to the accession.gbk.
my $input = setup_input( \%options );
my $count;

foreach my $id (@$input) {
    if ( $options{batch} ) {
        if ( $count < 500 ) {
            $count++;
        }
        if ( $count > 500 ) {
            sleep(8);
            $count = 0;
        }
    }
    my $efetch = Bio::DB::EUtilities->new(
        -db      => 'nucleotide',
        -id      => $id,
        -rettype => $rettype,
        -retmode => 'text',
        -email   => 'mymail@foo.bar',
    );

    my $out_prefix = defined $options{output_prefix} ? $options{output_prefix} : $id;
    open OUT, ">$out_dir/$out_prefix$out_suffix" or die "Unable to open $out_dir$out_prefix$out_suffix\n";
    print OUT $efetch->get_Response->content;
    close OUT;
    print STDERR "Downloaded: $id to: $out_dir/$out_prefix$out_suffix\n";
    if ( $rettype eq "gb" )    { &check_gb_file("$out_dir/$out_prefix$out_suffix"); }
    if ( scalar(@$input) > 1 ) { sleep 3; }                                             # Required so NCBI does not get mad.
}

sub check_gb_file {
    my $file = shift;
    chomp( my $tail = `tail -n 3 $file | head -n 1 ` );
    if ( $tail =~ /^CONTIG.+join\(([A-Za-z0-9\.]+)\:/ ) {
        my $second_id = $1;
        return if ( -e "$out_dir/$second_id$out_suffix" );
        print STDERR "The genbank file looks \"incomplete\". Downloading $second_id also.\n";
        my $efetch = Bio::DB::EUtilities->new(
            -db      => 'nucleotide',
            -id      => $second_id,
            -rettype => $rettype,
            -retmode => 'text',
            -email   => 'mymail@foo.bar',
        );
        open OUT, ">$out_dir/$second_id$out_suffix" or die "Unable to open $out_dir$second_id$out_suffix\n";
        print OUT $efetch->get_Response->content;
        close OUT;
        print STDERR "Downloaded: $second_id to: $out_dir/$second_id$out_suffix\n";
    }
}
