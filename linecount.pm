package linecount;
use strict;
use warnings;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( wc );

sub wc {
   my $file = shift;
   my $count=`wc -l $file`;
   $count=~/^(\d+)\s+/;
   my $num=$1;
   if($num==0){
      print STDERR "Warning: The file is empty. File:$file\tLines:$num\n";
   }
   return "$num";
}
1;
