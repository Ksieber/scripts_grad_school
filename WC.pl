#!/usr/bin/perl
use strict;
use warnings;
use linecount;

my $linecount = wc($ARGV[0]);
print "$linecount\n";
