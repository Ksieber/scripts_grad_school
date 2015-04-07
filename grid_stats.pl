#!/usr/bin/perl
use warnings;
use strict;

=head1 NAME

grid_stats.pl

=head1 SYNOPSIS

Generate a list of users and jobs they are running on the grid.

=head1 DESCRIPTION

Generate a list of users and jobs they are running on the grid.
May also pass @ARGV SGE user names to ignore.

=head1 EXAMPLE

watch -n 120 grid_stats.pl UserNameToQuery

=head1 AUTHOR - Karsten Sieber

e-mail: Karsten.sieber@gmail.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with "_"

=cut

if ( defined $ARGV[0] and ( $ARGV[0] =~ /\-{0,2}[?hH]{1}/ or $ARGV[0] =~ /help/i ) ) {
    die "\n\tThis script will list the number of jobs and threads for each SGE user. Any user name given via ARGV is ignored.\n\tExample: grid_stats.PL Gridhog1\n\n";
}

my $users_to_ignore = "";
if ( defined @ARGV ) {
    foreach my $user (@ARGV) {
        $users_to_ignore = $users_to_ignore . " | grep -v \"$user\"";
    }
}

my $stats = {};

my @qstat_lines = `qstat -u "*"$users_to_ignore`;
foreach my $line (@qstat_lines) {
    $line =~ /\d+\.\d+\s+\w+\s+(\w+)\s+(\w+)/;
    next if !$1;
    my $user       = $1;
    my $state      = $2;
    my @split_line = split( /\s+/, $line );
    my $threads    = ( $split_line[-1] =~ /\d+\-\d+\:\d+/ or $split_line[-1] =~ /\d{3,}/ ) ? $split_line[-2] : $split_line[-1];
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
