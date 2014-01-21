#!/usr/bin/perl
use strict;
use warnings;
use linecount;

my $linecount = wc($ARGV[0]);
if($linecount==0){
	print "This file is empty: $linecount.\n";
}
print "$linecount\n";
