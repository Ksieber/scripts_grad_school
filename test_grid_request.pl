#!/usr/bin/perl
use warnings;
use strict;
use Time::SoFar;
use Grid::Request;

my @array = ( 'Superman', 'Batman', 'Spiderman', 'X-men' );
chomp(@array);
my $request = Grid::Request->new(
    project    => "jdhotopp-lab",
    initialdir => "/home/ksieber/",
    error      => "/home/ksieber/",
    opsys      => "Linux",
);

$request->command("/local/projects-t3/HLGT/ksieber/tmp/echo.sh");
$request->add_param(
    {   type  => "ARRAY",
        key   => '$(Name)',
        value => \@array,
    }
);

$request->submit_and_wait();

#my $elapsed = runtime();
#print STDERR "$elapsed\n";
