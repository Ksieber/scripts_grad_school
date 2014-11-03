#!/usr/bin/perl
use warnings;
use strict;

my $num = 0;
while (<>) {
	$num=$num+$_;
}

print "$num\n";