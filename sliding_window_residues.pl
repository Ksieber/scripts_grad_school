#!/usr/local/bin/perl 
use strict;
use Bio::SeqIO;
use Getopt::Long;
use POSIX;

## Prints left most position of the window

my $infile = "";
my $winsize = 25;  #default
my $stepsize = 1; #default
GetOptions('input=s' => \$infile, 'winsize=i' => \$winsize, 'stepsize=i' => \$stepsize);
die("Usage = resisudes.pl --input=<fasta file> [--winsize=<int>] [--stepsize=<int>]\n") if(!$infile);

my $seqio = Bio::SeqIO->new('-file' => $infile, '-format' => 'fasta');
while(my $seqobj = $seqio->next_seq) {
    my $id  = $seqobj->display_id;
    my $len = $seqobj->length;

	my $grandtot;
	for(my $i = 1; $i <= ($len-($winsize-1)); $i+=$stepsize) {
		my $seq = $seqobj->subseq($i,$i+($winsize-1));
		my $tot = $seq =~ tr/a-zA-Z\*//;
		$grandtot += $tot;
		my $t = $seq =~ tr/tT//;
		my $c = $seq =~ tr/cC//;
		my $g = $seq =~ tr/gG//;
		my $a = $seq =~ tr/aA//;
		my $n = $seq =~ tr/aAcCgGtT//c;
		my $p = (($c+$g)/$tot)*100;
		my $purine = (($a+$g)/$tot)*100;
		print "$id\t$i\tperG+C:$p\tperA+G:$purine\tT:$t\tC:$c\tG:$g\tA:$a\n";
	}
}