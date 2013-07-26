package Range;
use strict;
use warnings;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( Range );


sub Range {
   shift if UNIVERSAL::isa( $_[ 0 ], __PACKAGE__ );
   my $arrayref = shift;
   return ( undef, undef ) unless defined $arrayref && @$arrayref > 0;
   my ( $zzz, $hi, $lo );
   $hi = $lo = $$arrayref[ 0 ];
   foreach $zzz ( @$arrayref ) {
      if ( $zzz > $hi ) {
         $hi = $zzz;
      }
      if ( $zzz < $lo ) {
         $lo = $zzz;
      }  
   }
   if ( $lo eq "" ) { $lo = "0" }
   return ( $hi, $lo );
}

1;
