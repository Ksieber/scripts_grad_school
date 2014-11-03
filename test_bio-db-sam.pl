#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
use warnings;
use strict;
use Carp;
use Bio::DB::Fasta;
use Bio::DB::Sam;
use Bio::Graphics;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use File::Basename;
use run_cmd;
use print_call;
use bwa;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
	'input|i=s',
	'ref=s',
	'region=s',
	'help|h'
) or confess "Invalid arguement. Please try agian.\n";

if($options{help}){&help;}
if(!$options{ref} || !$options{input} || !$options{region}) { confess "Error: You must pass all of the folowing arguements: --ref --bam --region. Please try agian.\n"; }
print_call(\%options);

my ($fn,$path,$suff) = fileparse($options{input},qr/\.[^\.]+/);
$options{region} =~ /(.*)\:(\d+)\-(\d+)/;
my ($chr,$start,$end) = ($1,$2,$3);

## Open input bams, creating index if needed
my $bam = Bio::DB::Sam->new(
	-bam => $options{input},
	-fasta => $options{ref1},
	-expand_flags => 1,
	-autoindex => 1);

my @aln = $bam->get_features_by_location(
	-seq_id => $chr,
	-start => $start,
	-end   => $end
	);

#print STDERR Dumper(@aln);

for my $a (@aln){
	my $seqid = $a->seq_id;
	my $start_aln = $a->start;
	my $end_aln   = $a->end;
	my $ref_dna = $a->dna;

	#print STDERR Dumper($seqid,$start_aln,$end_aln,$ref_dna);
}

my @pairs = $bam->get_features_by_location(
	-type => 'read_pair',
	-seq_id => $chr,
	-start => $start,
	-end   => $end
	);

#print STDERR Dumper(@pairs);
for my $pair (@pairs) {
	my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
	# print STDERR Dumper($first_mate,$second_mate);
	my $f_start = $first_mate->start;
	my $s_start = $second_mate->start; 
	# print STDERR Dumper($f_start,$s_start);
}


sub help {
	die "Help: This script will take a 2 bams mapped @ different references and merge the reference & map at the merged reference. This is useful for merging LGT regions.
		--input=				bam1. 5' side of merger. Assumes position sorted.
		--region=			<chr#:100-200>
		--ref1=				Reference 1 fasta.
		--help|h\n";
}
__END__



