package errorcheck;
use strict;
use warnings;
use Carp;
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw( errchk );

sub errchk {
    my $error = $_[0];
    if ( $error != 0 ) {
        confess "Error\n";
    }
}
1;
