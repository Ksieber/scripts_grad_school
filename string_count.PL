#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
use strict;
use warnings;

if ( !@ARGV ) { &help; }
if ( $ARGV[0] eq '--help' or $ARGV[0] eq '-h' or $ARGV[0] eq '-?' ) { &help; }

my $string = $ARGV[0];
print STDOUT length($string) . "\n";

sub help { die "This script will count the character length of the string passed via $ARGV[0]\n"; }
