#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::EUtilities;
my @ids = qw(CAB02640 EAS10332 YP_250808 NP_623143 P41007);

my $factory = Bio::DB::EUtilities->new(
    -eutil   => 'efetch',
    -db      => 'protein',
    -id      => \@ids,
    -email   => 'mymail@foo.bar',
    -rettype => 'gi'
);

my @gis = split( m{\n}, $factory->get_Response->content );

print join( ',', @gis ), "\n";
