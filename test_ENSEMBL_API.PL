#!/user/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $regfeat_adaptor = $registry->get_adaptor( 'Human', 'Funcgen', 'RegulatoryFeature' );

my $slice_adaptor   = $registry->get_adaptor( 'Human', 'Core',    'Slice' );
my $slice = $slice_adaptor->fetch_by_region( 'chromosome', 6, 133396267, 133397174 );

my @reg_feats = @{ $regfeat_adaptor->fetch_all_by_Slice($slice) };
foreach my $rf (@reg_feats) {
    print $rf->stable_id . ": ";
    print_feature($rf);
    print "\tCell: " . $rf->cell_type->name . "\n";
    print "\tFeature Type: " . $rf->feature_type->name . "\n";
}





print STDOUT "Pass\n";

#Prints absolute coordinates and not relative to the slice
sub print_feature {
    my $feature = shift;
    print $feature->display_label . "\t(" . $feature->seq_region_name . ":" . $feature->seq_region_start . "-" . $feature->seq_region_end . ")\n";
}


__END__
use HTTP::Tiny;
use JSON;

my $ext = $in_lower < $in_upper ? "/map/human/GRCh37/$in_chr\:$in_lower\..$in_upper\:1/GRCh38\?" : "/map/human/GRCh37/$in_chr\:$in_lower\..$in_upper\:-1/GRCh38\?";
my $response = $http->get( $server . $ext, { headers => { 'Content-type' => 'application/json' } } );
die "Failed!\n" unless $response->{success};

if ( length $response->{content} ) {
	my $hash = decode_json( $response->{content} );
            foreach my $mappings ( @{ $hash->{mappings} } ) { ... }
