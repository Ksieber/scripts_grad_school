#!/usr/bin/perl
use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,
    'input|i=s',    #
    'output_dir|o=s',
    'output_prefix=s',
    'evalue_cutoff=s',
    'krona|k=i',
    'best_hits_only|B=i',
    'Qsub|q=i',
    'help|h',
) or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) { &help; }

use lib "/home/ksieber/scripts";
use run_cmd;
use print_call;
use MultipleFileHandles;

if ( $options{Qsub} ) { &Qsub_script( \%options ) }

use Carp;
$Carp::MaxArgLen = 0;
use File::Basename;
use Bio::DB::EUtilities;
use Bio::DB::Taxonomy;
my $db = Bio::DB::Taxonomy->new( -source => 'entrez' );

print_call( \%options );
my ( $fn, $dir, $suf ) = fileparse( $options{input}, qr/\.[^\.]+/ );
my $out_dir = $options{output_dir} ? $options{output_dir} : $dir;
run_cmd("mkdir -p $out_dir");
my $out_pref = $options{output_prefix} ? $options{output_prefix} : "$fn\_lca";
my $lca_fh;
if ( $options{output_dir} || $options{output_prefix} ) {
    open( $lca_fh, ">", "$out_dir\/$out_pref\.txt" ) or confess "Error: Can't open output: $out_dir\/$out_pref\.txt because: $!\n";
}
else {
    $lca_fh = *STDOUT;
}
my $OUTPUT = MultipleFileHandles->new( [$lca_fh] );
if ( $options{krona} ) {
    open( my $krona_fh, "| cut -f2 | perl -ne 's/;/\\t/g; print;' | /home/ksieber/bin/ktImportText - -q -o $out_dir\/$out_pref\_krona.html" )
        or confess "Error: Can't open krona output: $out_dir\/$out_pref\_krona.html because: $!";
    $OUTPUT->add($krona_fh);
}
open( IN, "<", "$options{input}" ) or die "Can't open: $options{input} because: $!\n";

my $last_id;
my $lca;
my $bho           = $options{best_hits_only} ? $options{best_hits_only} : undef;
my $evalue_cutoff = $options{evalue_cutoff}  ? $options{evalue_cutoff}  : 1;
while (<IN>) {
    chomp;
    my @fields     = split( /\t/, $_ );
    my $current_id = $fields[0];
    my $evalue     = $fields[19];
    $fields[15] =~ /\[([A-Za-z0-9\s]+)\]/;
    my $organism = $1;
    next if ( $organism !~ /\w+/ );
    if ( !$last_id ) {
        $last_id = $current_id;
        $lca     = &organism2lineage($organism);
        if ($bho) { $evalue_cutoff = $evalue; }
    }
    if ( $last_id eq $current_id && $evalue <= $evalue_cutoff ) {
        my $lineage = &organism2lineage($organism);
        my $new_lca = &find_lca( [ $lca, $lineage ] );
        $lca = $new_lca;
    }
    if ( $current_id ne $last_id ) {
        $OUTPUT->print("$last_id\t$lca\n");
        $last_id = $current_id;
        $lca     = &organism2lineage($organism);
        if ($bho) { $evalue_cutoff = $evalue; }
    }
}

# Print the final LCA out.
$OUTPUT->print("$last_id\t$lca\n");

$OUTPUT->close;

## Complete!
print_complete( \%options );

sub help {
    die "Help: This script will take a btab output file and create the LCA for each query.
	--input|i=				btab blast.txt
	--output_dir|o=			Directory for output.										[*STDOUT]
	  --output_prefix=		Prefix_lca.txt	  											[Input filename]
	--krona|k=				<0|1> 1= Draw Krona plot. 									[0]
	--evalue_cutoff=		Max evalue.													[1]
	--best_hits_only|B=		<0|1> 1= Parse the btab blast report for best hits only. 	[0] 
	--Qsub|q=				<0|1> 1= qsub the job to the grid. 							[0]
	  --project=			Grid project to use. 										[jdhotopp-lab]
	--help|h

	Example: 	
		perl /home/ksieber/scripts/btab2lca.pl --input=input_blast.txt --output_dir=/some/where/for/output/ --krona=1 --evalue_cutoff=1e-10
		perl /home/ksieber/scripts/btab2lca.pl -i input_blast.txt -o /some/where/for/output/ -k 1 --best_hits_only=1
				\n";
}

sub organism2lineage {
    my $organism = shift;
    my $taxon    = $db->get_taxon( -name => $organism ) or return -1;
    my $name     = $taxon->scientific_name;
    my $c        = $taxon;
    my @lineage  = ($name);
    while ( my $parent = $db->ancestor($c) ) {
        unshift @lineage, $parent->scientific_name;
        $c = $parent;
    }
    my $lineage = join( ";", @lineage );
    return $lineage;
}

sub find_lca {
    my $lineages = shift;

    # prime LCA
    my @lca = split( ';', $lineages->[0] );

    foreach my $l (@$lineages) {
        my $newlca = [];
        my @lineage = split( ';', $l );
        for ( my $i = 0; $i < @lineage; $i++ ) {
            if ( $lca[$i] eq $lineage[$i] ) {
                push( @$newlca, $lineage[$i] );
            }
            else {
                last;
            }
        }
        @lca = @$newlca;
    }
    return join( ';', @lca );
}

__END__


 





