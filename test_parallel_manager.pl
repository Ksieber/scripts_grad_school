#!/usr/bin/perl
use warnings;
use strict;
use Time::SoFar;
use Parallel::ForkManager;


my $pm = Parallel::ForkManager->new(4);
my @superheroes = ('Superman','Batman','Spiderman','X-men');

foreach my $hero (@superheroes){
	$pm -> start and next;
	sleep 10;
	print STDERR "$hero\n";
	$pm -> finish;
}
$pm->wait_all_children;
my $elapsed = runtime();
print STDERR "Script completed in: $elapsed\n";