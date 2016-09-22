#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
use strict;
use warnings;

use Bio::DB::EUtilities;
use Data::Dumper;

my @ids = qw(257153399);

my $eutil = Bio::DB::EUtilities->new(
    -eutil => 'esummary',
    -email => 'mymail@foo.bar',
    -db    => 'nucleotide',
    -id    => \@ids
);

# $eutil->print_all;  ## Use this to find the _get_Items_by_name that you may want (ie the field "Title" in this)

while ( my $docsum = $eutil->next_DocSum ) {
    my ($item) = $docsum->get_Items_by_name('Title');
    print $item->get_content."\n";
}

