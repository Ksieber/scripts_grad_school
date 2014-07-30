#!/usr/bin/perl
use strict;
use warnings;
use read_in_list;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
	'opt1=s',
	'opt2=s',
	);

my $hash = {
	'foo' => '1',
	'bar' => '2',
	'salad' => '100',
}; 

print STDERR &find_key_with_min_hash_value($hash)."\n";

sub find_key_with_min_hash_value {
	my $hash   = shift;
    my ($key, @keys) = keys   %$hash;
    my ($small, @vals) = values %$hash;

    for (0 .. $#keys) {
        if ($vals[$_] < $small) {
            $small = $vals[$_];
            $key = $keys[$_];
        }
    }
    $key
}
