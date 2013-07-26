#!/usr/local/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
      'fasta=s',
      'fastq',
      'help',
      );
if(!$options{fasta} && !$ARGV[0]){die
   "Error. Please give an input fasta file to add \\1 & \\2 to the ID's. Use either --fasta=<input> or pass the input on command line (perl script.pl <input>).Can also use --fastq to turn on fastq format.\n";
}

my $format = $options{fastq} ? "fastq" : "fasta" ;
my $input = $options{fasta} ? $options{fasta} : $ARGV[0];
my %seen;

my $reader = new Bio::SeqIO(-format=>$format,-file=>$input);
while(my $entry=$reader->next_seq){
   my $header = $entry->display_id();
   my $seq = $entry->seq();
   $seen{$header}++;
   if($seen{$header} == 1){
      print "\>$header\/1\n$seq\n";
   } elsif ($seen{$header} == 2){
      print "\>$header\/2\n$seq\n";
      delete($seen{$header});
   }
}
__END__

