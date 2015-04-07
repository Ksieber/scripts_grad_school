#!/usr/bin/perl
use strict;
use warnings;

if ( -t STDIN and not @ARGV) {
    die "\nThis script will convert .fq to .fa. Please pass input to ARGV[0] or through STDIN.\n\n";
}
my $IN;
if ( defined $ARGV[0] and -e $ARGV[0] ) {
    open( $IN, "<", $ARGV[0] ) or die "Error: Unable to open input fromt ARGV[0]: $ARGV[0]\n";
}
else {
    $IN = *STDIN;
}

while (<$IN>) {
    if ( $_ !~ /\w+/ ) { die "Error: Unable to read input.\n"; }
    chomp( my $id  = $_ );
    chomp( my $seq = <$IN> );
    my $comment = <$IN>;
    my $qual    = <$IN>;
    $id =~ s/^\@/\>/;
    print STDOUT "$id\n$seq\n";
}
