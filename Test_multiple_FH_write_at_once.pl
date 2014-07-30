#!/usr/bin/perl
## Test printing to multiple FH @ once.
use warnings;
use strict;
use lib '/home/ksieber/scripts/';
use MultipleFileHandles;

open(my $foo, "> /local/projects-t3/HLGT/TCGA/ksieber_dir/tmp/test_mfh/Test_multi.tmp") or confess "This didn't work.";
open(my $bar, "> /local/projects-t3/HLGT/TCGA/ksieber_dir/tmp/test_mfh/Test_multi_other.tmp") or confess "This didn't work.";

my $mfh = MultipleFileHandles->new([$foo,$bar]);

open(my $fh3, "> foobarsalad.txt");
$mfh->add($fh3);

$mfh->print("THIS WORKED!\n");

$mfh->close;
