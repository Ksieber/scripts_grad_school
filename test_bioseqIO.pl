#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $gb_file = "/local/projects-t3/HLGT/TCGA/ksieber_dir/COAD/SRP009542/jdh_work/references/NC_003454-Fusobacterium_nucleatum.gb";
my $stream = Bio::SeqIO->new(-file => $gb_file,
							-format => 'GenBank');

while(my $seq = $stream->next_seq()) {
	my @feats = $seq->get_all_SeqFeatures();
	foreach my $feat (@feats) {
		my $primary = $feat->primary_tag();
		if($feat->has_tag('product')) {
			my @vals = $feat->get_tag_values('product');
			my @locus = $feat->get_tag_values('locus_tag');
			my $start = $feat->start;
			my $stop = $feat->end;
			my $obj = {
				'primary' => $primary,
				'product' => $vals[0],
				'locus' => $locus[0],
				'start' => $feat->start,
				'end' => $feat->end
                };
            print STDERR "$obj->{locus}\n";
        }
    }
}