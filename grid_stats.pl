#!/usr/bin/perl

use warnings;
use strict;

my $stats = {};

my @qstat_lines = `qstat -u "*"`;
foreach my $line (@qstat_lines) {
    $line =~ /\d+\.\d+\s+\w+\s+(\w+)\s+(\w+)/;
    next if !$1;
    my $user       = $1;
    my $state      = $2;
    my @split_line = split( /\s+/, $line );
    my $threads    = $split_line[-1];
    $stats->{$user}->{total}++;
    $stats->{$user}->{state}->{$state}->{$threads}++;
}

foreach my $user ( sort { $stats->{$b}->{total} <=> $stats->{$a}->{total} } keys %$stats ) {
    printf "%-15s %-10s", $user, $stats->{$user}->{total};
    foreach my $states ( sort { $b cmp $a } keys %{ $stats->{$user}->{state} } ) {
        foreach my $threads ( sort keys %{ $stats->{$user}->{state}->{$states} } ) {
            printf " %-10s", "$stats->{$user}->{state}->{$states}->{$threads}\:$states\:$threads";
        }
    }
    print "\n";
}
