 #!/usr/bin/perl -ws

use Bio::DB::Taxonomy;
my $db = Bio::DB::Taxonomy->new(-source => 'entrez');




while ( <> ) {
	my $line = $_;
	chomp $line;
	# tab-del file with <organism>	<count>
	my @f = split (/\t/, $line);
	
	# use NCBI Entrez over HTTP
	my $taxonid = $db->get_taxonid($f[0]);

	# get a taxon
	my $taxon = $db->get_taxon(-taxonid => $taxonid);

	if(!$taxon) {
		print STDERR "Couldn't find taxonomy info for $f[0]\n";
		$retval = {
			'acc'=> $acc,
			'gi'=> $gi,
			'taxon_id'=> $taxonid
		};
	}
	elsif ($taxon->isa('Bio::Taxon')) {
		my $name = $taxon->scientific_name;
		my $c = $taxon;
		my @lineage = ($name);
		while (my $parent = $db->ancestor($c)) {
			unshift @lineage, $parent->scientific_name;
			$c = $parent;
		}
		$retval = {
			'gi' => $gi,
			'acc'=> $acc,
			'taxon_id'=> $taxonid, 
			'name'=> $name, 
			'lineage'=> join(";", @lineage)
		};
		print $f[1], "\t", $f[0], "\t", $retval->{lineage}, "\t", $retval->{gi}, "\n";
	}else {
		print STDERR "Had something other than a Bio::Taxon\n";
	}
}
