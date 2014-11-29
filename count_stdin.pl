#!/usr/bin/perl
use warnings;
use strict;

my $num = 0;
while (<>) {
    if ( $_ =~ /\w+\s{1}\:\s{1}(\d+)$/ ) {
        $num += $1;
    }
    else {
        $num += $_;
    }
}

print "$num\n";
