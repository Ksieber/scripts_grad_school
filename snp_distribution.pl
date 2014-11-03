#!/usr/bin/perl

use strict;
use warnings;
use Carp;
$Carp::MaxArgLen = 0;
use POSIX;

my $file = $ARGV[0];
if(!$file){confess "Error: Must pass a vcf file to snp_distribution.pl\n";}

my %contig_length;
my %snp_distribution;

open(IN, "<","$file") or confess "Can't open input: $file\n";
while(<IN>){
	chomp;
	my $line=$_;
	if($line=~/^\#\#contig/){
		$line =~/^##.*\<ID=(.*),length=(\d+)\>$/;
		my $contig = $1;
		my $length = $2;
		$contig_length{$contig}=$length;
	} elsif ($line=~/^\#/){
		next;
	} else {
		my ($contig,$position,$info_col,$ten_col)=(split/\t/,$line)[0,1,7,9];
		next if ($info_col=~/INDEL/);
		my @f=split(/:/,$ten_col);
		next if ($f[0]!~/1\/1/ || $f[2]<99);
		my $percent = ceil(($position/$contig_length{$contig})*100);
		$snp_distribution{$percent}++;
	}
}

for(my $i=1; $i<=100; $i++){
	print "$i\t$snp_distribution{$i}\n";
}
