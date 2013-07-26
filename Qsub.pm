package Qsub;
use strict;
use warnings;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( Qsub );


sub Qsub {
   my $shell_script = shift;
   my $report=`qsub -V -P jdhotopp-lab $shell_script`;
   print STDERR "$report";
}
1;
