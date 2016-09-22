#!/usr/bin/perl
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use strict;
use warnings;
use Carp;
use mk_dir;
use setup_input;
use read_in_list;
use File::Basename;
use Statistics::R;
use Data::Dumper;

if ( !@ARGV ) { &help_full; }

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,    'input|i=s',   'input_list|I=s', 'abundance|a=s',     'top_n_otus=i', 'column|c=i',      'lca_cut=i',       'key=i',
    'na_color=s', 'output|O=s',  'output_dir|o=s', 'output_prefix|p=s', 'cluster|C=i',  'cluster_num|n=i', 'bottom_margin=i', 'right_margin=i',
    'help|h',     'help_full|?', 'linker_file|L=s',
) or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} )      { &help; }         ## &help is @ the end of the script
if ( $options{help_full} ) { &help_full; }    ## &help is @ the end of the script

if ( !$options{input} and !$options{input_list} ) { die "Error: Must pass an input of LCA.txt to convert to a matrix using --input= or --input_list=[files.list] or [file1.txt,file2.txt]\n"; }

my $inputs = setup_input( \%options );
my ( $fn, $path, $dir ) = defined $options{input} ? fileparse( $options{input}, qr/\.[^\.]+/ ) : fileparse( $options{input_list}, qr/\.[^\.]+/ );
my $output_dir    = defined $options{output_dir}    ? $options{output_dir}    : $path;
my $output_prefix = defined $options{output_prefix} ? $options{output_prefix} : $fn;
my $output        = defined $options{output}        ? $options{output}        : "$output_dir/$output_prefix\_heatmap.pdf";
$output =~ s/\/{2,}/\//g;                     ## Cleanup output name with multiple "/" in it. ie: /path//for//output.pdf -> /path/for/output.pdf
my $column        = defined $options{column}        ? $options{column}        : 1;
my $key           = defined $options{key}           ? $options{key}           : "1";
my $top_n_otus    = defined $options{top_n_otus}    ? $options{top_n_otus}    : "20";
my $bottom_margin = defined $options{bottom_margin} ? $options{bottom_margin} : "16";
my $right_margin  = defined $options{right_margin}  ? $options{right_margin}  : "12";
my $lca_cut       = defined $options{lca_cut}       ? $options{lca_cut}       : "7";
my $na_color      = defined $options{na_color}      ? $options{na_color}      : "#999999";
my $cluster       = defined $options{cluster}       ? $options{cluster}       : 0;
if ( !$options{abundance} ) { $options{abundance} = "LN"; }
$options{abundance} =~ tr/[a-z]/[A-Z]/;

my $links;
if ( $options{linker_file} ) {
    open( my $LINKER, "<", $options{linker_file} ) or die "Error: Unable to open the linker_file: $options{linker_file} because: $options{linker_file}\n";
    while (<$LINKER>) {
        chomp( my $line = $_ );
        next if ( $line =~ /^#/ );
        my @split_line = split( /\s+/, $line );
        $links->{ $split_line[0] } = $split_line[1];
    }
    close $LINKER;
}

my $data;
my $total_otu_per_id;
my %unique_otu_counts;
my $otu_ranks;

# 1st, read in the LCA's from all the files
my %id_hash;
foreach my $input (@$inputs) {
    my ( $fn, $path, $suf ) = fileparse( $input, ( '_bwa-lca_independent_lca.txt', '_bwa-lca.txt', '_lca-bwa.txt', 'blastn-m8_lca-independent.out', qr/\.[^\.]+/ ) );
    open( my $IN, "< $input" ) or die "Error: Unable to open the input: $input\n";
    while (<$IN>) {
        chomp;
        my @split_line = ( split /\t/, $_ );
        next if ( scalar(@split_line) == 1 );
        my $read_id   = $split_line[0];
        my $lca       = $split_line[$column];
        my @split_lca = split( /;/, $lca );
        my @short_lca = ( $lca_cut eq "-1" ) ? @split_lca : splice( @split_lca, 0, $lca_cut );
        my $otu       = $short_lca[-1];
        $otu =~ s/\s+/__/g;
        $otu =~ s/\//_/g;
        $otu =~ s/\-/_/g;
        $otu =~ s/\(//g;
        $otu =~ s/\)//g;
        $otu =~ s/\'//g;
        $otu =~ s/\"//g;
        $otu =~ s/\.//g;

        if ( defined $options{linker_file} ) {
            if ( !$otu ) { die "Foo1:$_"; }
            if ( !$links->{$read_id} ) { die "Foo2:$read_id\n"; }
            $data->{ $links->{$read_id} }->{$otu}++;
            $total_otu_per_id->{ $links->{$read_id} }++;
            $id_hash{ $links->{$read_id} }++;
        }
        else {
            $data->{$fn}->{$otu}++;
            $total_otu_per_id->{$fn}++;
            $id_hash{$fn}++;
        }
        $unique_otu_counts{$otu}++;
        $otu_ranks->{$otu} = scalar(@short_lca);    # Wrough estimate of specificity of the OTU (higher rank = more specific)
    }
    close $IN;

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
my $num_of_files = scalar( keys %id_hash );
my $rownames_string = "\"" . join( "\", \"", ( sort { $a cmp $b } keys %id_hash ) ) . "\"";

# 5th, Print some nice text for the user
print STDERR "\n\t*** lca2heatmap.pl ***\n";
if ( $options{abundance} eq 'RELATIVE' or $options{abundance} eq 'REL' ) { print STDERR "\t*** Creating heatmap of the relative abundance. ***\n"; }
elsif ( $options{abundance} eq 'LOG' )    { print STDERR "\t*** Creating heatmap of the LOG transformed relative abundance. ***\n"; }
elsif ( $options{abundance} eq 'LN' )     { print STDERR "\t*** Creating heatmap of the LN transformed relative abundance. ***\n"; }
elsif ( $options{abundance} eq 'COUNTS' ) { print STDERR "\t*** Creating heatmap of the counts of OTU abundance. ***\n"; }
else                                      { die "*** ERROR *** Unable to determine the abundance transform: $options{abundance}\n"; }

# 6th, load & manipulate the data into R.
my $R = Statistics::R->new( r_bin => '/home/ksieber/bin/R', tmp_dir => '/local/projects-t3/HLGT/TCGA/ksieber_dir/tmp/' );
$R->run('library(gplots)');
$R->run("reds = rev(rainbow(200))[1:10]");
$R->run("COLORS = c(rainbow(200)[32:200],reds)");
my $first_data = 1;

foreach my $otu (@top_n_otus_list) {
    my @list_of_otus_per_id;

    foreach my $id ( sort { $a cmp $b } keys %id_hash ) {
        if ( defined $data->{$id}->{$otu} ) {
            if ( $options{abundance} eq 'LOG' ) {
                my $log = log( ( $data->{$id}->{$otu} / $total_otu_per_id->{$id} ) * 100 ) / log(10);
                push( @list_of_otus_per_id, $log );
            }
            if ( $options{abundance} eq 'LN' ) {
                my $ln = log( ( $data->{$id}->{$otu} / $total_otu_per_id->{$id} ) * 100 );
                push( @list_of_otus_per_id, $ln );
            }
            elsif ( $options{abundance} eq 'RELATIVE' or $options{abundance} eq 'REL' ) {
                my $relative = ( $data->{$id}->{$otu} / $total_otu_per_id->{$id} ) * 100;
                push( @list_of_otus_per_id, $relative );
            }
            elsif ( $options{abundance} eq 'COUNTS' ) {
                push( @list_of_otus_per_id, $data->{$id}->{$otu} );
            }
        }
        else {
            push( @list_of_otus_per_id, "NA" );
        }
    }

    my $otus_per_id_string = join( ", ", @list_of_otus_per_id );
    $R->run("$otu<-matrix(c($otus_per_id_string), nrow=$num_of_files, ncol=1)");

    # Cluster doesn't allow NA, but we want NA in there to color by later. Make a seperate matrix to cluster by with 0 instead of NA
    if ( $cluster == 1 ) {
        my $otus_per_id_string_for_clustering = $otus_per_id_string;
        $otus_per_id_string_for_clustering =~ s/(NA)/0/g;
        $R->run("$otu\_for_clustering<-matrix(c($otus_per_id_string_for_clustering), nrow=$num_of_files, ncol=1)");
        if   ( $first_data == 1 ) { $R->run("data_for_cluster<-cbind($otu\_for_clustering)"); }
        else                      { $R->run("data_for_cluster<-cbind(data_for_cluster,$otu\_for_clustering)"); }
    }

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

# 7th, cluster the data (if $options{cluster}=1)
my $R_string_dendrogram;
my $R_string_row_clustering;
my $R_string_col_clustering;
if ( $cluster == 1 ) {
    my $max_number_of_clusters = scalar( keys %id_hash );
    my $number_of_clusters = ( defined $options{cluster_num} ) ? $options{cluster_num} : "5";
    if ( $number_of_clusters > $max_number_of_clusters ) {
        print STDERR " *** Warning: You have too many clusters. Try lowering --clusters_num=. Currently setting to $num_of_files.\n";
        $number_of_clusters = $max_number_of_clusters;
    }
    ## Order the data by col Sums
    $R->run("DataOrd = data_for_cluster[,order(colSums(data_for_cluster),decreasing=\"T\")]");
    ## Cluster the data
    $R->run("DataClustered = hclust(dist(as.matrix(DataOrd)), method=\"ward.D\")");

    ## Break the cluster into number of desired clusters
    $R->run("DataTree = cutree(DataClustered,k=$number_of_clusters)");
    $R_string_dendrogram     = "row";
    $R_string_row_clustering = "as.dendrogram(DataClustered)";
    $R_string_col_clustering = "colnames(data)";
}
else {
    $R_string_dendrogram     = "none";
    $R_string_row_clustering = "NA";
    $R_string_col_clustering = "colnames(data)";
}

# 8th, Actually graph the LCAs
$R->run("pdf(file=\"$output\")");
$R->run(
    "heatmap.2(data, 
    col=COLORS,
    na.color=\"$na_color\",
    dendrogram=\"$R_string_dendrogram\",
    Rowv=$R_string_row_clustering,
    Colv=$R_string_col_clustering,
    trace=\"none\",
    margin=c($bottom_margin,$right_margin),
    keysize=1,
    key=$key,
    )"
);
$R->run("dev.off()");
$R->stopR();

print STDERR "\t*** OUTPUT = $output ****\n\n";

sub help_full {
    die "This script will create a HEATMAP from file(s) with: \"{read_name}    {LCA}\" file.
    --input|i=          LCA File.
    --input_list|I=     List of LCA files, either a file with 1 file/line or a comma seperated list.
      --column=         < # >   [1]  Column to select LCA from the input, 0 based counting. 
    --linker_file|L=    < /path/to/linker.txt > Tab delimited: read-id and unique-id (participant, analysis, etc). 
    --output|O=         < /path/and/name/for/output.pdf >
      --output_dir|o=   < /path/for/output/ >
      --output_prefix|p={prefix_for_output}.pdf
    --abundance=        < LOG | LN | RELATIVE | COUNTS > [LN] Transform abundance. 
                LOG & LN transform relative abundance. REL shorthand for RELATIVE. 
    --cluster|C=        <0|1>   [0] 1= Cluster similier samples.
    --cluster_num=      < # >   [5] Number of clusters to create.
    --top_n_otus=       < # >   [20] Draw the \${top_n_otus}
    --lca_cut=          < # >   [7]  Cut the LCA Taxonomy at this split #. (These #s only work for ~microbiome).
       0=Life 1=Domain 2=Kingdom 3=Phylum 4=Class 5=Order 6=Family 7=Genus 8=Species 9=Strain -1=All
    --na_color=         [\"#999999\"] Color to substitute for NA data. Takes any \"color\" R:heatmap2 can accept.
    --key=              < 0|1 > [1]  Draw Heatmap key.
    --bottom_margin=    < # >   [16] Heatmap bottom_margin.
    --right_margin=     < # >   [12] Heatmap right_margin.
    --help|h 
    --help_full|?\n";
}
__END__

