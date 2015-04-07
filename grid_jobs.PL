#!/usr/bin/perl
use warnings;
use strict;

=head1 NAME

grid_jobs.PL

=head1 SYNOPSIS

Generate a list of the number of each SGE job-name for a user(s)

=head1 EXAMPLE

watch -n 120 grid_jobs.PL ksieber GridHog1

=head1 AUTHOR - Karsten Sieber

e-mail: Karsten.Sieber@gmail.com

=cut

if ( !@ARGV ) {
    my $user = `whoami`;
    chomp($user);
    push( @ARGV, $user );
}

if ( $ARGV[0] =~ /\-{0,2}[?hH]{1}/ or $ARGV[0] =~ /help/i ) {
    die "\n\tThis script will make a list of the number of SGE job-names running per user passed via ARGV. Default=user running script.\n\tExample: grid_jobs.PL ksieber Gridhog1\n\n";
}

foreach my $users (@ARGV) {
    my @qstat_lines = `qstat -u "$users"`;
    my %number_of_jobs;
    foreach my $line (@qstat_lines) {
        chomp($line);

        next if ( $line =~ /^job-ID/ );
        next if ( $line =~ /^[-]{2,}/ );

        $line =~ /\d+\s+[0-9\.]+\s+([A-Za-z0-9\-\_]+)/;
        my $job_name = $1;
        $number_of_jobs{$job_name}++;

    }
    print "\# $users jobs:\n";
    foreach my $job_name ( sort { $number_of_jobs{$b} <=> $number_of_jobs{$a} } keys %number_of_jobs ) {
        printf "\t%-10s %-5s", $number_of_jobs{$job_name}, $job_name;
        print "\n";
    }
}

