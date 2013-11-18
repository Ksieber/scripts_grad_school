#!/usr/local/bin/perl 
use strict;
use Bio::SeqIO;
use Getopt::Long;

my $infile = "";
my $winsize = 7;  #default
GetOptions('infile=s' => \$infile, 'winsize=i' => \$winsize, 'stepsize=i' => \$stepsize);
die("Usage = window.pl --infile  <fasta file> [--winsize <int>] [--stepsize <int>]\n") if(!$infile);
$stepsize = defined $stepsize ? $stepsize : 1;

my $seqio = Bio::SeqIO->new('-file' => $infile, '-format' => fasta);
while(my $seqobj = $seqio->next_seq) {
    my $id  = $seqobj->display_id;
    my $len = $seqobj->length;

for(my $i = 1; $i <= $len-($winsize-1)); $i+=$stepsize) {
    my $window = $seqobj->subseq($i,$i+($winsize-1));
    # do clever stuff with window
	}
}