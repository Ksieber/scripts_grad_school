#!/usr/bin/perl
use lib ( '/home/ksieber/scripts/' );
use strict;
use linecount;

if ( scalar(@ARGV) > 1 ) {
    foreach my $file (@ARGV) {
        if ( !-e $file ) { warn "This file doens't exist, skipping: $file\n"; }
        my $linecount = wc($file);
        print "$file : $linecount\n";
    }
}
elsif ( -e $ARGV[0] ) {
    my $linecount = wc( $ARGV[0] );
    print "$linecount\n";
}
elsif ( !@ARGV and !-t STDIN ) {
    my $counter;
    while (<STDIN>) {
        $counter++;
    }
    print "$counter\n";
}
