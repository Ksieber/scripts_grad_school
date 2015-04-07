package print_progress;
use warnings;
use strict;
our @ISA    = qw(Exporter);
our @EXPORT = qw( print_progress );
$| = 1;
my $counter;

sub print_progress {
    $counter++;
    if ( $counter >= 1000000 ) { $counter = 0; print STDERR "\r                                                                                        "; }
    elsif ( $counter % 100000 == 0 ) {
        my $string;
        for ( my $n = 1; $n <= ( $counter / 100000 ); $n++ ) { $string = $string . ".\t"; }
        print STDERR "\r$string";
    }
}

1;
