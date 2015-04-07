#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
use warnings;
use strict;
use Data::Dumper;
use XML::Simple qw(:strict);

my $xml = $ARGV[0];
my $ref = XMLin( $xml, KeyAttr => { Result => 'id' }, ForceArray => [ 'Result', 'file' ], GroupTags => { files => 'file' } )
    or die "======== &downloadCGHub: &XMLin is unable to read in the cgquery.xml file: $xml because: $!\n";

foreach my $sample ( sort keys %{ $ref->{'Result'} } ) {
    print STDERR "SAMPLE: " . Dumper($sample) . "\n";
    foreach my $file ( @{ $ref->{'Result'}->{$sample}->{files} } ) {
        print STDERR "FILE: " . Dumper($file) . "\n";
        print STDERR "FILENAME: $file->{filename}\tMD5: $file->{checksum}->{content}\n";
    }
}
