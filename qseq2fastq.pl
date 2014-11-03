#!/usr/bin/perl

use warnings;
use strict;

while (<>) {
      chomp;
      my @parts = split /\t/;
      next if($parts[10]==0);
      print "@","$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]#$parts[6]/$parts[7]\n";
      $parts[8] =~ tr/\./N/;
      print "$parts[8]\n";
      print "+\n";
#$parts[9] =~ tr/\x40-\xff\x00-\x3f/\x21-\xe0\x21/;
      print "$parts[9]\n";
}
