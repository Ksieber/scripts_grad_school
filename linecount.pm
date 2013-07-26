package linecount;
use strict;
use warnings;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( wc );

sub wc {
   my $file = $_[0];
   my $count=`wc -l $file`;
   $count=~/^(\d+)\s+/;
   my $num=$1;
   if($num!~/\d+/){
      print STDERR "Warning: The file didn't have any lines?\n";
   }
   return "$num";
}
1;
