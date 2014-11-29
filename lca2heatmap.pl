#!/usr/bin/perl
use strict;
use warnings;
use lib ( '/home/ksieber/perl5/lib/perl5/', '/home/ksieber/scripts/' );
use Carp;
use mk_dir;
use setup_input;
use read_in_list;
use File::Basename;
use Statistics::R;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results
    = GetOptions( \%options, 'input|i=s', 'input_list|I=s', 'top_n_otus=i', 'column|c=i', 'lca_cut=i', 'key=i', 'output_dir|o=s', 'output_prefix|p=s', 'bottom_margin=i', 'right_margin=i', 'help|?' )
    or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} ) { &help; }    ## &help is @ the end of the script

if ( !$options{input} and !$options{input_list} ) { die "Error: Must pass an input of LCA.txt to convert to a matrix using --input= or --input_list=[files.list] or [file1.txt,file2.txt]\n"; }

my $inputs = setup_input( \%options );
my ( $fn, $path, $dir ) = defined $options{input} ? fileparse( $options{input}, qr/\.[^\.]+/ ) : fileparse( $options{input_list}, qr/\.[^\.]+/ );
my $output_dir    = defined $options{output_dir}    ? $options{output_dir}    : $path;
my $output_prefix = defined $options{output_prefix} ? $options{output_prefix} : $fn;
my $output        = "$output_dir/$output_prefix\_heatmap.pdf";
my $column        = defined $options{column}        ? $options{column}        : 1;
my $key           = defined $options{key}           ? $options{key}           : "1";
my $top_n_otus    = defined $options{top_n_otus}    ? $options{top_n_otus}    : "20";
my $bottom_margin = defined $options{bottom_margin} ? $options{bottom_margin} : "16";
my $right_margin  = defined $options{right_margin}  ? $options{right_margin}  : "12";
my $lca_cut       = defined $options{lca_cut}       ? $options{lca_cut}       : "7";
my $data;
my %unique_otu_counts;
my $otu_ranks;

# 1st, read in the LCA's from all the files
my @filename_list;
foreach my $input (@$inputs) {
    my ( $fn, $path, $suf ) = fileparse( $input, ( '_bwa-lca_independent_lca.txt', '_bwa-lca.txt', '_lca-bwa.txt', 'blastn-m8_lca-independent.out', qr/\.[^\.]+/ ) );
    open( my $IN, "< $input" ) or die "Error: Unable to open the input: $input\n";
    while (<$IN>) {
        chomp;
        my @split_line = ( split /\t/, $_ );
        next if ( scalar(@split_line) == 1 );
        my $lca       = $split_line[$column];
        my @split_lca = split( /;/, $lca );
        my @short_lca = splice( @split_lca, 0, $lca_cut );
        my $otu       = $short_lca[-1];
        $otu =~ s/\s+/__/g;
        $otu =~ s/\//_/g;
        $otu =~ s/\-/_/g;
        $otu =~ s/\(//g;
        $otu =~ s/\)//g;
        $otu =~ s/\'//g;
        $otu =~ s/\"//g;
        $otu =~ s/\.//g;
        $data->{$fn}->{$otu}++;
        $unique_otu_counts{$otu}++;
        $otu_ranks->{$otu} = scalar(@split_lca);    # Wrough estimate of specificity of the OTU (higher rank = more specific)
    }
    close $IN;
    push( @filename_list, $fn );
}

# 2nd, generate a list of the top X number of OTUs and create a string for the colnames in R.
my @top_n_otus_list;
my $colname_string;
if ( $top_n_otus > 0 ) {

    # Sort high to low the number of times each OTU was seen
    my @high_to_low_otu_counts = sort { $b <=> $a } values %unique_otu_counts;
    my $otu_count_cutoff = ( scalar( keys %unique_otu_counts ) < ( $top_n_otus - 1 ) ) ? $high_to_low_otu_counts[-1] : $high_to_low_otu_counts[ $top_n_otus - 1 ];

    # Grab the top N otu's seen OTU key
    my %top_otu_keys;

    foreach my $otu ( keys %unique_otu_counts ) {
        if ( $unique_otu_counts{$otu} >= $otu_count_cutoff ) {
            $top_otu_keys{$otu}++;
        }
    }

    # Sort the list of top N OTUs based on their rank score
    @top_n_otus_list = sort { $otu_ranks->{$a} <=> $otu_ranks->{$b} } keys %top_otu_keys;
    $colname_string = "\"" . join( "\", \"", @top_n_otus_list ) . "\"";
}
else {
    $colname_string = "\"" . join( "\", \"", ( keys %$otu_ranks ) ) . "\"";
}

# 4th, calculate the number of filenames we have and create a string for the rownames in R.
my $num_of_files = scalar(@filename_list);
my $rownames_string = "\"" . join( "\", \"", @filename_list ) . "\"";

# 5th, load & manipulate the data in R.
my $R = Statistics::R->new( r_bin => '/home/ksieber/bin/R' );
$R->run('library(gplots);');
# $R->run('library(heatmap3)');
my $first_data = 1;

foreach my $otu (@top_n_otus_list) {
    my @otus_per_fn_list;
    foreach my $fn (@filename_list) {
        if ( defined $data->{$fn}->{$otu} ) {
            push( @otus_per_fn_list, $data->{$fn}->{$otu} );
        }
        else {
            push( @otus_per_fn_list, "0" );
        }
    }
    my $otus_per_fn_string = join( ", ", @otus_per_fn_list );
    $R->run("$otu<-matrix(c($otus_per_fn_string), nrow=$num_of_files, ncol=1)");
    if ( $first_data == 1 ) {
        $R->run("data<-cbind($otu)");
        $first_data++;
    }
    else {
        $R->run("data<-cbind(data,$otu)");
    }
}
$R->run("rownames(data)<-c($rownames_string)");
$R->run("colnames(data)<-c($colname_string)");

$R->run("reds = rev(rainbow(200))[1:10]");
$R->run("COLORS = c(rainbow(200)[32:200],reds)");

$R->run("pdf(file=\"$output\")");
$R->run(
    "heatmap.2(data, 
	col=COLORS,
	dendrogram=\"none\",
	Rowv=NA,
	Colv=NA,
	trace=\"none\",
	margin=c($bottom_margin,$right_margin),
	keysize=1,
	key=$key,
	)"
);
$R->run("dev.off()");

__END__

