#!/usr/local/bin/perl -w
use strict;
use average;
use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
      'fasta=s',
      
      'help',
      );
if(!$options{fasta} && !$ARGV[0]){die
   "Error. Please give an input file to calculate the min and max length on. Use either --fasta=<input> or pass the input on command line (perl script.pl <input>).\n";
}

my $max;
my $min;
my $total_bp;
my $total_reads;

my $input = $options{fasta} ? $options{fasta} : $ARGV[0];

my $reader = new Bio::SeqIO(-format=>'fasta',-file=>$input);
while(my $seqRec=$reader->next_seq){
   my $length = $seqRec->length;
   if(!$max){$max=$length};
   if(!$min){$min=$length};
   if($length > $max){$max=$length;}
   if($length < $min){$min=$length;}   
   $total_bp += $length;
   $total_reads++;
}

my $average = $total_bp / $total_reads;

print "Done processing $input.\nMax length= $max\nMin length= $min\nAverage: $average\n";

__END__

