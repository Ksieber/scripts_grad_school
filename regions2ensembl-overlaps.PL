#!/usr/bin/perl
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use strict;
use warnings;
use setup_input;
use Bio::EnsEMBL::Registry;
use Cwd;
use mk_dir;
use File::Basename;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'input_list|I=s', 'ensembl|e=s', 'output|O=s', 'output_prefix|p=s', 'output_dir|o', 'help|?' )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( defined $options{help} or not @ARGV ) { &help; }
if ( !defined $options{input} and !defined $options{input_list} ) {
    die "Error: Must pass a region to query using --input=<chr:start-end> or a list of regions with --input_list=</path/to/file.list>\n";
}

my @ensembl_query_split = defined $options{ensembl} ? split( /,/, $options{ensembl} ) : ( 'regs', 'transcripts', 'genes', 'repeats' );
my %ensembl_query_hash;
foreach my $queries (@ensembl_query_split) { $ensembl_query_hash{$queries}++; }

my $regs        = defined $ensembl_query_hash{'regs'}        ? 1 : 0;
my $transcripts = defined $ensembl_query_hash{'transcripts'} ? 1 : 0;
my $genes       = defined $ensembl_query_hash{'genes'}       ? 1 : 0;
my $repeats     = defined $ensembl_query_hash{'repeats'}     ? 1 : 0;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);
my $slice_adaptor   = $registry->get_adaptor( 'Human', 'Core',    'Slice' );
my $regfeat_adaptor = $registry->get_adaptor( 'Human', 'Funcgen', 'RegulatoryFeature' );
my $tr_adaptor      = $registry->get_adaptor( 'Human', 'Core',    'Transcript' );

my %inputs;
if ( defined $options{input} ) {
    $inputs{ $options{input} }++;
}
if ( defined $options{input_list} ) {
    open( IN, "<", $options{input_list} ) or die "Error: Unable to open in input_list:$options{input_list}\n";
    my $i = 0;
    while (<IN>) { next if ( $_ =~ /^\#/ ); $i++; my @split_line = split( /\s+/, $_ ); $inputs{ $split_line[0] } = $i; }
    close IN;
}

my ( $fn, $dir, $suf ) = defined $options{input_list} ? fileparse( $options{input_list}, qr/\.[^\.]+/ ) : ( "tmp_fn", getcwd, ".txt" );
my $out_file = defined $options{output} ? $options{output} : "$dir$fn\_regions2ensembl.txt";
my ( $out_fn, $out_dir, $out_suf ) = fileparse( $out_file, qr/\.[^\.]+/ );
mk_dir($out_dir);
my $OUT_FH;

if ( defined $options{output} or defined $options{output_dir} or defined $options{output_prefix} ) {
    open( $OUT_FH, ">", "$out_file" ) or die "Error: Unable to open the output: $out_file\n";
}
else {
    $OUT_FH = *STDOUT;
}

foreach my $region ( sort { $inputs{$a} <=> $inputs{$b} } keys %inputs ) {
    $region =~ /^([A-Za-z0-9\|\.]+)\:(\d+)\-(\d+)/;
    my ( $chr, $lower, $upper ) = ( $1, $2, $3 );
    if    ( $chr =~ /chr([0-9XY]+)/ ) { $chr = $1; }
    elsif ( $chr =~ /^([0-9XY]+)$/ )  { $chr = $1; }
    else                              { next; }

    printf $OUT_FH ( "%-30s\t", $region );
    my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $lower, $upper );

    # Genes
    if ( $genes == 1 ) {
        print $OUT_FH "genes:";
        print_features( $slice->get_all_Genes() );
    }

    # Transcripts
    if ( $transcripts == 1 ) {
        print $OUT_FH "transcripts:";
        print_features( $tr_adaptor->fetch_all_by_Slice($slice) );
    }

    # Regulatory Features
    if ( $regs == 1 ) {
        print $OUT_FH "regs:";
        my @hit_slices = @{ $regfeat_adaptor->fetch_all_by_Slice($slice) };
        my $rf_hit     = 0;
        my %non_redudant_hits;

        my $i = 0;
        while ( my $rf = shift @hit_slices ) {
            $i++;
            $non_redudant_hits{ $rf->feature_type->name } = $i;
        }

        if ( scalar( keys %non_redudant_hits ) >= 1 ) {
            print $OUT_FH join( ";", ( sort { $non_redudant_hits{$a} <=> $non_redudant_hits{$b} } keys %non_redudant_hits ) ) . "\t";
        }
        else { print $OUT_FH "Na\t"; }
    }

    # Repeats
    if ( $repeats == 1 ) {
        print $OUT_FH "repeats:";

        my @repeats    = @{ $slice->get_all_RepeatFeatures() };
        my $repeat_hit = 0;
        my %non_redudant_hits;

        my $i = 0;
        while ( my $repeat = shift @repeats ) {
            $i++;
            $non_redudant_hits{ $repeat->display_id() } = $i;
        }

        if ( scalar( keys %non_redudant_hits ) >= 1 ) {
            print $OUT_FH join( ";", ( sort { $non_redudant_hits{$a} <=> $non_redudant_hits{$b} } keys %non_redudant_hits ) ) . "\t";
        }
        else { print $OUT_FH "Na\t"; }
    }

    print $OUT_FH "\n";
}

sub print_features {
    my $features = shift;
    my %non_redudant_hits;

    my $i = 0;
    while ( my $hit = shift @{$features} ) {
        $i++;
        $non_redudant_hits{ $hit->display_xref->display_id() } = $i;
    }
    if ( scalar( keys %non_redudant_hits ) >= 1 ) {
        print $OUT_FH join( ";", ( sort { $non_redudant_hits{$a} <=> $non_redudant_hits{$b} } keys %non_redudant_hits ) ) . "\t";
    }
    else {
        print $OUT_FH "Na\t";
    }
}

sub help {
    die "
    This script will pull info from ENSEMBL based on regions from hg38. 
    --input|i=              <chr:start-end>
    --input_list|I=         </path/to/file/regions.list
    --ensembl|e=            Default is all these ('regs', 'transcripts', 'genes', 'repeats'). User may specify >=1 comma seperated.

    Output goes to STDOUT unless a --output* option is used.
    --output|O=             </full/output/path.txt>
    --output_dir|o=         </dir/for/outputs/>
    --output_prefix=        <{prefix}.txt>\n\n";
}
