package average;
use strict;
use warnings;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( Mean );


sub Mean {
   shift if UNIVERSAL::isa( $_[ 0 ], __PACKAGE__ );
   my $arrayref = shift;
   return undef unless defined $arrayref && @$arrayref > 0;
   my $result;
   foreach ( @$arrayref ) { $result += $_ }
   return $result / @$arrayref;
}
1;
