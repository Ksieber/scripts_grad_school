#!/usr/bin/perl
use strict;
use warnings;

if ( !defined $ARGV[0] ) { die "Error: This script must have ARGV0 to set the new shell TMPDIR.\n"; }
print "export TMPDIR=$ARGV[0]\n";