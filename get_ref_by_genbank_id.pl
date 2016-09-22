#!/usr/bin/perl -I /home/ksieber/scripts/ -I /home/ksieber/perl5/lib/perl5/

=head1 NAME

get_ref_by_genbank_id.pl

=head1 SYNOPSIS

USAGE:  get_ref_by_genbank_id.pl --input=NC_0000001

=head2 OPTIONS

=over 24

=item B<--input|acc|i=>

Single accession

=item B<--input_list|list|I=>

List of accessions, either comma seperate or a file with 1 per line.

=item B<--type=>

"FASTA", "gb", or "both"; Default is FASTA.

=item B<--output=>

/full/path/and/name/for/output.fa

=item B<--output_dir|o=>

/dir/path/for/output/

=item B<--output_prefix|p=>

${output_prefix}.fa

=back

=head1 DESCRIPTION

get_ref_by_genbank_id.pl - Retrieves genbank records from genbank.

=cut
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use strict;
use warnings;
use Bio::DB::EUtilities;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Basename;
use setup_input;
use Cwd;
use mk_dir;
use print_call;

# Process the options list.
my %options;
my $results = GetOptions( \%options, 'input|acc|i=s', 'input_list|list|I=s', 'type|t=s', 'batch', 'help|h|?', 'output|O=s', 'output_dir|o=s', 'output_prefix|p=s', ) or &pod2usage();

&pod2usage( -exitval => 1, -verbose => 2, -output => \*STDOUT ) if ( $options{help} or ( !$options{input} and !$options{input_list} ) );

my @rettype = ( defined $options{type} and $options{type} =~ m/both/i ) ? ( "FASTA", "gb" ) : ( defined $options{type} ? $options{type} : 'FASTA' );
my $out_dir = defined $options{output_dir} ? $options{output_dir} : getcwd;
$out_dir =~ s/\/{2,}/\//g;
$out_dir =~ s/\/$//;
mk_dir($out_dir);

# Loop through the IDs and pull down the genbank file.
# Write the results to the accession.gbk.
my $input = setup_input( \%options );
my $count;

print STDERR "\n";

foreach my $id (@$input) {
    foreach my $rettype (@rettype) {
        $rettype =~ s/fasta/FASTA/;    ## force caps for EUtils
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
        my $out_suffix = ( $rettype eq 'FASTA' )         ? ".fa"                   : ".gbk";
        my $retry      = 1;
    GET_SEQ: my $response = $efetch->get_Response->content;
        if ( $response =~ /Resource temporarily unavailable/ ) {
            if ( $retry == 1 ) {
                $retry = 0;
                sleep 90;
                goto GET_SEQ;
            }
            else {
                die "*** FAILED to download: $id ***\n";
            }
        }
        open OUT, ">$out_dir/$out_prefix$out_suffix" or die "Unable to open $out_dir$out_prefix$out_suffix\n";
        print OUT $response;
        close OUT;
        print STDERR "\tDownloaded: $id to: $out_dir/$out_prefix$out_suffix\n";
        if ( $rettype eq "gb" )    { &check_gb_file("$out_dir/$out_prefix$out_suffix"); }
        if ( scalar(@$input) > 1 ) { sleep 30; }                                            # Required so NCBI does not get mad.
    }
}

print STDERR "\n";
print_complete( \%options );

sub check_gb_file {
    my $file = shift;
    my ( $out_prefix, $out_dir, $out_suffix ) = fileparse( $file, ".gbk" );
    chomp( my $tail = `tail -n 3 $file | head -n 1 ` );
    if ( $tail =~ /^CONTIG.+join\(([A-Za-z0-9\.\_]+)\:/ ) {
        my $second_id = $1;
        return if ( -e "$out_dir/$second_id$out_suffix" );
        print STDERR "\n\tThe genbank file may be \"incomplete\". Downloading $second_id also.\n\n";
        my $efetch = Bio::DB::EUtilities->new(
            -db      => 'nucleotide',
            -id      => $second_id,
            -rettype => "gb",
            -retmode => 'text',
            -email   => 'mymail@foo.bar',
        );
        open OUT, ">$out_dir/$second_id$out_suffix" or die "Unable to open $out_dir/$second_id$out_suffix\n";
        print OUT $efetch->get_Response->content;
        close OUT;
        print STDERR "\tDownloaded: $second_id to: $out_dir/$second_id$out_suffix\n";

        if ( defined $options{type} and $options{type} =~ m/both/i ) {
            my $efetch_fasta = Bio::DB::EUtilities->new(
                -db      => 'nucleotide',
                -id      => $second_id,
                -rettype => "FASTA",
                -retmode => 'text',
                -email   => 'mymail@foo.bar',
            );
            open OUT, ">$out_dir/$second_id\.fa" or die "Unable to open $out_dir/$second_id\.fa\n";
            print OUT $efetch_fasta->get_Response->content;
            close OUT;
            print STDERR "\tDownloaded: $second_id to: $out_dir/$second_id\.fa\n";
        }
    }
}
