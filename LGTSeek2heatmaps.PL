#!/usr/local/bin/perl
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use strict;
use warnings;
use Carp;
use mk_dir;
use Cwd;
use File::HomeDir;
use POSIX;
use print_call;
use run_cmd;
use setup_input;
use read_in_list;
use File::Basename;
use Data::Dumper;
use Statistics::R;
if ( !@ARGV ) { &help_full; }
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,       'microbiome_lca|m=s', 'lgt_lca|l=s',     'abundance|a=s',     'top_n_otus=i',  'top_otus_by=s',   'column|c=i',         'lca_cut=i',
    'key=i',         'na_color=s',         'output_dir|o=s',  'output_prefix|p=s', 'cluster|C=i',   'cluster_num|n=s', 'bottom_margin=i',    'right_margin=i',
    'help|h',        'help_full|?',        'linker_file|L=s', 'cluster_by=s',      'metadata|M=s@', 'ColorTable=i',    'SidebarColorSize=i', 'ColorTable_help',
    'cluster_den=i', 'sort_bac=i',         'sort_bac_by=s',   'nboot=i',           'threads|t=i',   'input=s',         'Qsub|q=i',           'sub_mail=i',
    'sub_mem=s',     'sub_name=s'
) or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} )            { &help; }               ## @ the end of the script
if ( $options{help_full} )       { &help_full; }          ## @ the end of the script
if ( $options{ColorTable_help} ) { &ColorTable_help; }    ## @ the end of the script
&verify_inputs;                                           ## @ the end of the script

my $out_dir = defined $options{output_dir} ? $options{output_dir} : getcwd;
mk_dir($out_dir);
my $out_prefix   = defined $options{output_prefix} ? $options{output_prefix} : "LGTSeek2heatmaps";
my $MICRO_output = "$out_dir/$out_prefix\_Microbiome-heatmap.pdf";
my $LGT_output   = "$out_dir/$out_prefix\_LGT-heatmap.pdf";
$MICRO_output =~ s/\/{2,}/\//g;
$LGT_output =~ s/\/{2,}/\//g;
my $column        = defined $options{column}        ? $options{column}        : 1;
my $key           = defined $options{key}           ? $options{key}           : 1;
my $top_n_otus    = defined $options{top_n_otus}    ? $options{top_n_otus}    : 20;
my $top_otus_by   = defined $options{top_otus_by}   ? $options{top_otus_by}   : "BOTH";
my $bottom_margin = defined $options{bottom_margin} ? $options{bottom_margin} : 16;
my $right_margin  = defined $options{right_margin}  ? $options{right_margin}  : 12;
my $lca_cut       = defined $options{lca_cut}       ? $options{lca_cut}       : 7;
my $na_color      = defined $options{na_color}      ? $options{na_color}      : "#999999";
my $cluster       = defined $options{cluster}       ? $options{cluster}       : 0;
if ( defined $options{cluster_by} ) { $cluster = 1; }
my $cluster_den = defined $options{cluster_den} ? $options{cluster_den} : 1;
my $nboot       = defined $options{nboot}       ? $options{nboot}       : 1000;
my $sort_bac    = defined $options{sort_bac}    ? $options{sort_bac}    : 0;
if ( ( defined $options{sort_bac_by} ) and ( !defined $options{sort_bac} ) ) { $sort_bac = 1; }
my $sort_bac_by      = defined $options{sort_bac_by}      ? $options{sort_bac_by}      : "TOTAL";
my $SidebarColorSize = defined $options{SidebarColorSize} ? $options{SidebarColorSize} : 2;
my $ColorTable       = defined $options{ColorTable}       ? $options{ColorTable}       : 1;
my $total_threads    = defined $options{threads}          ? $options{threads}          : 3;
$options{threads} = $total_threads;
my $threads     = $total_threads - 1;
my $cluster_fxn = "hclust";
if ( !defined $options{abundance} ) { $options{abundance} = "LN"; }
$options{abundance} =~ tr/[a-z]/[A-Z]/;
if ( !defined $options{cluster_by} ) { $options{cluster_by} = "MICRO"; }
$options{cluster_by} =~ tr/[a-z]/[A-Z]/;

if ( defined $options{Qsub} and $options{Qsub} == 1 ) {
    $options{input}      = $options{lgt_lca};
    $options{output_dir} = $out_dir;
    if ( !defined $options{sub_name} ) { $options{sub_name} = "LGTSeek2Heatmap"; }
    if ( !defined $options{sub_mem} )  { $options{sub_mem}  = ceil( $threads * .75 ) . "G"; }
    $options{metadata} = join( ',', @{ $options{metadata} } );
    Qsub_script( \%options );
}

print_call( \%options );
print STDERR "\n\t*** LGTSeek2heatmaps.PL ***\n";

# Load the Linker data
print STDERR "\t*** Reading in Linker data: $options{linker_file} ***\n";
my $links;
open( my $LINKER, "<", $options{linker_file} ) or die "Error: Unable to open the linker_file: $options{linker_file} because: $options{linker_file}\n";
while (<$LINKER>) {
    chomp( my $line = $_ );
    next if ( $line =~ /^#/ );
    my @split_line = split( /\s+/, $line );
    $links->{ $split_line[0] } = $split_line[1];
}
close $LINKER;

# Load the LCA data
my $otu_ranks;
print STDERR "\t*** Reading in microbiome data: $options{microbiome_lca} ***\n";
my $MICRO_data = &read_in_lca_file( $options{microbiome_lca} );
print STDERR "\t*** Reading in lgt data: $options{lgt_lca} ***\n";
my $LGT_data = &read_in_lca_file( $options{lgt_lca} );
my $merged_data = &merged_micro_lgt_data( $MICRO_data, $LGT_data );

my @top_n_otus_list;
my $colname_string;
my %top_MICRO_otu_keys;
my %top_LGT_otu_keys;
my %top_BOTH_otu_keys;

print STDERR "\t*** Determining which lcas to illustrate ***\n";
if ( $top_n_otus > 0 ) {
    $top_otus_by =~ tr/[a-z]/[A-Z]/;
    if ( $top_otus_by eq 'TOTAL' ) {

        # In case the user specified --top_n_otus > # of otus we have data for
        my $otu_count_cutoff
            = ( scalar( keys %{ $merged_data->{'otus_total'} } ) < $top_n_otus ) ? scalar( keys %{ $merged_data->{'otus_total'} } ) : $top_n_otus;

        # Grab the top N otu's
        my %top_otu_keys;
        my $counter = 1;

        for my $otu ( sort { $merged_data->{'otus_total'}->{$b} <=> $merged_data->{'otus_total'}->{$a} } keys %{ $merged_data->{'otus_total'} } ) {
            last if ( $counter > $otu_count_cutoff );
            $top_otu_keys{$otu}++;
            $counter++;
        }

        # Sort the list of top N OTUs based on their rank score
        @top_n_otus_list = sort { $otu_ranks->{$a} <=> $otu_ranks->{$b} } keys %top_otu_keys;
        $colname_string = "\"" . join( "\", \"", @top_n_otus_list ) . "\"";
    }
    else {

        if ( $top_otus_by eq 'MICRO' or $top_otus_by eq 'BOTH' ) {

            # In case the user specified --top_n_otus > # of otus we have data for
            my $otu_count_cutoff
                = ( scalar( keys %{ $MICRO_data->{'otus_total'} } ) < $top_n_otus ) ? scalar( keys %{ $MICRO_data->{'otus_total'} } ) : $top_n_otus;

            # Grab the top N otu's
            my %top_otu_keys;
            my $counter = 1;

            for my $otu ( sort { $MICRO_data->{'otus_total'}->{$b} <=> $MICRO_data->{'otus_total'}->{$a} } keys %{ $MICRO_data->{'otus_total'} } ) {
                last if ( $counter > $otu_count_cutoff );
                $top_MICRO_otu_keys{$otu}++;
                $counter++;

            }
        }
        if ( $top_otus_by eq 'LGT' or $top_otus_by eq 'BOTH' ) {

            # In case the user specified --top_n_otus > # of otus we have data for
            my $otu_count_cutoff
                = ( scalar( keys %{ $LGT_data->{'otus_total'} } ) < $top_n_otus ) ? scalar( keys %{ $LGT_data->{'otus_total'} } ) : $top_n_otus;

            # Grab the top N otu's
            my %top_otu_keys;
            my $counter = 1;

            for my $otu ( sort { $LGT_data->{'otus_total'}->{$b} <=> $LGT_data->{'otus_total'}->{$a} } keys %{ $LGT_data->{'otus_total'} } ) {
                last if ( $counter > $otu_count_cutoff );
                $top_LGT_otu_keys{$otu}++;
                $counter++;
            }
        }

        # Sort the list of top N OTUs based on their rank score
        if ( $top_otus_by eq 'BOTH' ) {
            map {
                my %hash = %{$_};
                foreach my $otu ( keys %hash ) { $top_BOTH_otu_keys{$otu}++; }
            } ( \%top_MICRO_otu_keys, \%top_LGT_otu_keys );
            @top_n_otus_list = sort { $otu_ranks->{$a} <=> $otu_ranks->{$b} } keys %top_BOTH_otu_keys;
        }
        elsif ( $top_otus_by eq 'MICRO' ) {
            @top_n_otus_list = sort { $otu_ranks->{$a} <=> $otu_ranks->{$b} } keys %top_MICRO_otu_keys;
        }
        elsif ( $top_otus_by eq 'LGT' ) {
            @top_n_otus_list = sort { $otu_ranks->{$a} <=> $otu_ranks->{$b} } keys %top_LGT_otu_keys;
        }

        $colname_string = "\"" . join( "\", \"", @top_n_otus_list ) . "\"";
    }
}
else {
    @top_n_otus_list = keys %{ $merged_data->{otus_total} };
}
if ( $sort_bac == 1 ) {
    print STDERR "\t*** Sorting bac lcas ***\n";
    @top_n_otus_list = @{ &reorder_bac_by_counts( \@top_n_otus_list ); };
}
$colname_string = "\"" . join( "\", \"", @top_n_otus_list ) . "\"";

# 4th, calculate the number of filenames we have and create a string for the rownames in R.
my $num_of_ids = scalar( keys %{ $merged_data->{'ids'} } );
my $rownames_string = "\"" . join( "\", \"", ( sort { $a cmp $b } keys %{ $merged_data->{'ids'} } ) ) . "\"";

print STDERR "\t*** Calculating abundances and loading data into R ***\n";

# 6th, load & manipulate the data into R.
my $original_tmpdir = `echo \$TMPDIR`;
`eval \$(/home/ksieber/scripts/change_shell_tmpdir.pl /tmp)`;
$ENV{'TMPDIR'} = '/tmp';
my $R = Statistics::R->new( r_bin => '/home/ksieber/bin/R', tmp_dir => '/tmp' );
$R->run("library(snow)");

if ( $cluster_fxn eq "pvclust" or $cluster_fxn eq "both" ) { $R->run("cl<-makeCluster($threads\, type=\"MPI\", outfile=\"\")"); }
$R->run("library(pvclust)");
$R->run("source(\"/home/ksieber/scripts/dendroCol.R\")");
$R->run('source("/home/ksieber/scripts/heatmap.3.R")');
&initialize_colors_in_R;

my $first_data = 1;
foreach my $otu (@top_n_otus_list) {
    ## Load LGT data into R.
    my @list_of_MICRO_otus_per_id;
    my @list_of_LGT_otus_per_id;

    foreach my $id ( sort { $a cmp $b } keys %{ $merged_data->{'ids'} } ) {

        # Micro data First
        if   ( defined $MICRO_data->{'otu_by_id'}->{$id}->{$otu} ) { push( @list_of_MICRO_otus_per_id, $MICRO_data->{'otu_by_id'}->{$id}->{$otu} ); }
        else                                                       { push( @list_of_MICRO_otus_per_id, "NA" ); }

        # LGT data Second
        if   ( defined $LGT_data->{'otu_by_id'}->{$id}->{$otu} ) { push( @list_of_LGT_otus_per_id, $LGT_data->{'otu_by_id'}->{$id}->{$otu} ); }
        else                                                     { push( @list_of_LGT_otus_per_id, "NA" ); }
    }
    my $MICRO_otus_per_id_string = join( ", ", @list_of_MICRO_otus_per_id );
    my $LGT_otus_per_id_string   = join( ", ", @list_of_LGT_otus_per_id );

    $R->run("MICRO_$otu<-matrix(c($MICRO_otus_per_id_string), nrow=$num_of_ids, ncol=1)");
    $R->run("LGT_$otu<-matrix(c($LGT_otus_per_id_string), nrow=$num_of_ids, ncol=1)");

    if ( $first_data == 1 ) {
        $R->run("MICRO_data<-cbind(MICRO_$otu)");
        $R->run("LGT_data<-cbind(LGT_$otu)");
        $first_data++;
    }
    else {
        $R->run("MICRO_data<-cbind(MICRO_data,MICRO_$otu)");
        $R->run("LGT_data<-cbind(LGT_data,LGT_$otu)");
    }
}

$R->run("rownames(MICRO_data)<-c($rownames_string)");
$R->run("colnames(MICRO_data)<-c($colname_string)");
$R->run("rownames(LGT_data)<-c($rownames_string)");
$R->run("colnames(LGT_data)<-c($colname_string)");

# Transform Microbiome data

if ( $options{abundance} eq 'LOG' ) {
    $R->run("MICRO_data[is.na(MICRO_data)]<-0");
    $R->run("MICRO_4heatmap<-t(apply(MICRO_data, 1, function(x) log(100 * x/sum(x))/log(10)))");
    $R->run("MICRO_4heatmap[!is.finite(MICRO_4heatmap)]<-NA");
}
elsif ( $options{abundance} eq 'LN' ) {
    $R->run("MICRO_data[is.na(MICRO_data)]<-0");
    $R->run("MICRO_4heatmap<-t(apply(MICRO_data, 1, function(x) log(100 * x/sum(x))))");
    $R->run("MICRO_4heatmap[!is.finite(MICRO_4heatmap)]<-NA");
}
elsif ( $options{abundance} eq 'RELATIVE' or $options{abundance} eq 'REL' ) {
    $R->run("MICRO_data[is.na(MICRO_data)]<-0");
    $R->run("MICRO_4heatmap<-t(apply(MICRO_data, 1, function(x) 100 * x/sum(x)))");
}
elsif ( $options{abundance} eq 'COUNTS' ) { $R->run("MICRO_4heatmap<-MICRO_data"); }

# Transform LGT data

if ( $options{abundance} eq 'LOG' ) {
    $R->run("LGT_data[is.na(LGT_data)]<-0");
    $R->run("LGT_4heatmap<-t(apply(LGT_data, 1, function(x) log(100 * x/sum(x))/log(10)))");
    $R->run("LGT_4heatmap[!is.finite(LGT_4heatmap)]<-NA");
}
elsif ( $options{abundance} eq 'LN' ) {
    $R->run("LGT_data[is.na(LGT_data)]<-0");
    $R->run("LGT_4heatmap<-t(apply(LGT_data, 1, function(x) log(100 * x/sum(x))))");
    $R->run("LGT_4heatmap[!is.finite(LGT_4heatmap)]<-NA");
}
elsif ( $options{abundance} eq 'RELATIVE' or $options{abundance} eq 'REL' ) {
    $R->run("LGT_data[is.na(LGT_data)]<-0");
    $R->run("LGT_4heatmap<-t(apply(LGT_data, 1, function(x) 100 * x/sum(x)))");
}
elsif ( $options{abundance} eq 'COUNTS' ) { $R->run("LGT_4heatmap<-LGT_data"); }

# 7th, cluster the data (if $options{cluster}=1)
my $R_string_dendrogram;
my $R_string_row_clustering;
my $R_string_col_clustering = "colnames(MICRO_data)";
my $smart_cluster_num;

if ( $cluster == 1 ) {
    if    ( $options{cluster_by} eq 'MICRO' ) { $R->run("data_for_id_cluster<-MICRO_data"); }
    elsif ( $options{cluster_by} eq 'LGT' )   { $R->run("data_for_id_cluster<-LGT_data"); }

    # First, transform counts into proportion
    $R->run("data_for_id_cluster<-t(apply(data_for_id_cluster, 1, function(x) 100 * x/sum(x)))");
    $R->run("data_for_id_cluster[is.na(data_for_id_cluster)]<-0");

    # Cluster the proportion data
    if ( $cluster_fxn eq "hclust" ) {
        print STDERR "\t*** Clustering ids with hclust ***\n";
        $R->run("DataClustered = hclust(dist(as.matrix(data_for_id_cluster)), method=\"ward.D\")");
        if ( defined $options{cluster_num} and $options{cluster_num} eq "SMART" ) {
            print STDERR "\t*** Cluster# based on Calinski-Harabasz index ***\n";
            $R->run("require(clusterSim)");
            $R->run("nclusters=NULL");
            $R->run("for(n in 1:20){ if(n==1){ nclusters[n]=NA } else { nclusters[n]=index.G1(data_for_id_cluster, as.vector(cutree(DataClustered,k=n)), centrotypes=\"centroids\") }}");
            $R->run("opt_cluster_num <- which(nclusters==max(nclusters, na.rm=TRUE))");
            $R->run("DataTree = cutree(DataClustered,k=opt_cluster_num)");
            $smart_cluster_num = $R->get("opt_cluster_num");
        }
        elsif ( defined $options{cluster_num} and $options{cluster_num} =~ /\d+/ ) { $R->run("DataTree = cutree(DataClustered,k=$options{cluster_num})"); }
        else                                                                       { $R->run("DataTree = cutree(DataClustered,h=(max(DataClustered\$height)/1.5))"); }

        # Set global variables for drawing the heatmap
        $R_string_dendrogram     = "row";
        $R_string_row_clustering = "as.dendrogram(DataClustered)";
    }
    if ( $cluster_fxn eq "both" ) {
        print STDERR "\t*** Clustering ids with hclust ***\n";
        $R->run("Both_DataClustered = hclust(dist(as.matrix(data_for_id_cluster)), method=\"ward.D\")");
        if ( defined $options{cluster_num} and $options{cluster_num} eq "SMART" ) {
            print STDERR "\t*** Cluster# based on Calinski-Harabasz index ***\n";
            $R->run("require(clusterSim)");
            $R->run("nclusters=NULL");
            $R->run("for(n in 1:20){ if(n==1){ nclusters[n]=NA } else { nclusters[n]=index.G1(data_for_id_cluster, as.vector(cutree(Both_DataClustered,k=n)), centrotypes=\"centroids\") }}");
            $R->run("opt_cluster_num <- which(nclusters==max(nclusters, na.rm=TRUE))");
            $R->run("Both_DataTree = cutree(Both_DataClustered,k=opt_cluster_num)");
        }
        elsif ( defined $options{cluster_num} and $options{cluster_num} =~ /\d+/ ) { $R->run("Both_DataTree = cutree(Both_DataClustered,k=$options{cluster_num})"); }
        else                                                                       { $R->run("Both_DataTree = cutree(Both_DataClustered,h=(max(Both_DataClustered\$height)/1.5))"); }
    }
    if ( $cluster_fxn eq "pvclust" or $cluster_fxn eq "both" ) {
        print STDERR "\t*** Clustering ids with pvclust ***\n";

        $R->run("DataClustered <- parPvclust(cl, as.matrix(t(data_for_id_cluster)), nboot=$nboot)");
        if ( defined $options{cluster_num} and $options{cluster_num} eq "SMART" ) {
            print STDERR "\t*** Cluster# based on Calinski-Harabasz index ***\n";
            $R->run("require(clusterSim)");
            $R->run("nclusters=NULL");
            $R->run("for(n in 1:20){ if(n==1){ nclusters[n]=NA } else { nclusters[n]=index.G1(data_for_id_cluster, as.vector(cutree(DataClustered,k=n)), centrotypes=\"centroids\") }}");
            $R->run("opt_cluster_num <- which(nclusters==max(nclusters, na.rm=TRUE))");
            $R->run("DataTree = cutree(DataClustered,k=opt_cluster_num)");
            $smart_cluster_num = $R->get("opt_cluster_num");
        }
        elsif ( defined $options{cluster_num} and $options{cluster_num} =~ /\d+/ ) { $R->run("DataTree = cutree(DataClustered\$hclust,k=$options{cluster_num})"); }
        else                                                                       { $R->run("DataTree = cutree(DataClustered\$hclust,h=(max(DataClustered\$hclust\$height)/1.5))"); }

        # Color the dendrogram based on sig pvclustering
        $R->run("sig_cluster_keys <- unlist(pvpick(DataClustered, alpha=0.95, pv=\"au\", type=\"geq\", max.only=TRUE)\$clusters)");
        $R->run("dend_colored <- dendrapply(as.dendrogram(DataClustered\$hclust), dendroCol, keys=sig_cluster_keys, xPar=\"edgePar\", bgr=\"black\", fgr=\"red\", pch=20)");

        # Set global variables for drawing the heatmap
        $R_string_dendrogram     = "row";
        $R_string_row_clustering = "dend_colored";
    }
}
else {
    $R_string_dendrogram     = "none";
    $R_string_row_clustering = "NA";
}
if ( $cluster_den == 0 ) { $R_string_dendrogram = "none"; }
#########################################
## Set up sidebars
my $heatmap_first_bar = 1;
my @sidebar_names;
if ( defined $options{metadata} ) {
    foreach my $sidebar_files ( split( /,/, join( ',', @{ $options{metadata} } ) ) ) {
        if ( !-e $sidebar_files ) { die "Error: The id_metadata file doesn't exist: $sidebar_files\n"; }
        my $sidebar_first_data = 1;
        my $metadata_hash      = &hash_in_data($sidebar_files);
        foreach my $id ( sort { $a cmp $b } keys %{ $merged_data->{'ids'} } ) {
            if ( !$metadata_hash->{$id} ) { die "Error: This id: $id wasn't found in the metadata file: $sidebar_files.\n"; }
            if ( $sidebar_first_data == 1 ) { $R->run("SidebarData <- cbind($metadata_hash->{$id})"); $sidebar_first_data++; }
            else                            { $R->run("SidebarData <- cbind(SidebarData, $metadata_hash->{$id})"); }
        }
        if   ( $ColorTable == 1 ) { $R->run("SidebarData <- colorTbl[as.vector(SidebarData)]"); }
        else                      { $R->run("SidebarData <- as.vector(SidebarData)"); }

        push( @sidebar_names, $metadata_hash->{'#Header'} );

        if ( $heatmap_first_bar == 1 ) { $R->run("MergeAllSidebars <- cbind(SidebarData)"); $heatmap_first_bar++; }
        else                           { $R->run("MergeAllSidebars <- cbind(MergeAllSidebars,SidebarData)"); }
    }
}

if ( $cluster == 1 ) {
    if ( $cluster_fxn eq "both" ) {

        # Build the Clustering sidebar
        if ( ( !defined $options{cluster_num} ) or ( $options{cluster_num} > 15 ) ) {
            $R->run("big_Table<-sample(rainbow(100))");
            $R->run("Both_ClusteringSideBars <- big_Table[as.vector(Both_DataTree)]");
        }
        else { $R->run("Both_ClusteringSideBars <- colorTbl[as.vector(Both_DataTree)]"); }
        push( @sidebar_names, "hclust" );

        # Add the clustering sidebar to other potential sidebars
        if ( $heatmap_first_bar == 1 ) { $R->run("MergeAllSidebars <- cbind(Both_ClusteringSideBars)"); $heatmap_first_bar++; }
        else                           { $R->run("MergeAllSidebars <- cbind(MergeAllSidebars,Both_ClusteringSideBars)"); }
    }

    # Build the Clustering sidebar
    if ( ( !defined $options{cluster_num} ) or ( defined $options{cluster_num} and $options{cluster_num} =~ /\d+/ and $options{cluster_num} > 15 ) or ( $smart_cluster_num > 15 ) ) {
        $R->run("big_Table<-sample(rainbow(100))");
        $R->run("ClusteringSideBars <- big_Table[as.vector(DataTree)]");
    }
    else { $R->run("ClusteringSideBars <- colorTbl[as.vector(DataTree)]"); }

    if   ( ( $cluster_fxn eq "both" ) or ( $cluster_fxn eq "pvclust" ) ) { push( @sidebar_names, "pvclust" ); }
    else                                                                 { push( @sidebar_names, "Clustering" ); }

    # Add the clustering sidebar to other potential sidebars
    if ( $heatmap_first_bar == 1 ) { $R->run("MergeAllSidebars <- cbind(ClusteringSideBars)"); $heatmap_first_bar++; }
    else                           { $R->run("MergeAllSidebars <- cbind(MergeAllSidebars,ClusteringSideBars)"); }
}
my $sidebar_names_R_string = "\"" . join( "\", \"", @sidebar_names ) . "\"";
$R->run("colnames(MergeAllSidebars)=c($sidebar_names_R_string)");
$R->run("MergeAllSidebars<-as.matrix(t(MergeAllSidebars))");

# 5th, Print some nice text for the user
if ( $options{abundance} eq 'RELATIVE' or $options{abundance} eq 'REL' ) { print STDERR "\t*** Creating heatmap of the relative abundance. ***\n"; }
elsif ( $options{abundance} eq 'LOG' )    { print STDERR "\t*** Creating heatmap of the LOG transformed relative abundance. ***\n"; }
elsif ( $options{abundance} eq 'LN' )     { print STDERR "\t*** Creating heatmap of the LN transformed relative abundance. ***\n"; }
elsif ( $options{abundance} eq 'COUNTS' ) { print STDERR "\t*** Creating heatmap of the counts of OTU abundance. ***\n"; }
else                                      { die "*** ERROR *** Unable to determine the abundance transform: $options{abundance}\n"; }

if ( $cluster == 1 or defined $options{metadata} ) {

    # 8th, create the MICRO Heatmap
    $R->run("pdf(file=\"$MICRO_output\")");
    $R->run(
        "heatmap.3(MICRO_4heatmap, 
    col=COLORS,
    na.color=\"$na_color\",
    dendrogram=\"$R_string_dendrogram\",
    Rowv=$R_string_row_clustering,
    RowSideColors=MergeAllSidebars,
    RowSideColorsSize=$SidebarColorSize,
    Colv=$R_string_col_clustering,
    trace=\"none\",
    margin=c($bottom_margin,$right_margin),
    keysize=1,
    key=$key,
    )"
    );
    $R->run("dev.off()");

    # 9th, create the LGT Heatmap
    $R->run("pdf(file=\"$LGT_output\")");
    $R->run(
        "heatmap.3(LGT_4heatmap, 
    col=COLORS,
    na.color=\"$na_color\",
    dendrogram=\"$R_string_dendrogram\",
    Rowv=$R_string_row_clustering,
    RowSideColors=MergeAllSidebars,
    RowSideColorsSize=$SidebarColorSize,
    Colv=$R_string_col_clustering,
    trace=\"none\",
    margin=c($bottom_margin,$right_margin),
    keysize=1,
    key=$key,
    )"
    );
    $R->run("dev.off()");
}

else {
    # 8th, create the MICRO Heatmap
    $R->run("pdf(file=\"$MICRO_output\")");
    $R->run(
        "heatmap.3(MICRO_4heatmap, 
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

    # 9th, create the LGT Heatmap
    $R->run("pdf(file=\"$LGT_output\")");
    $R->run(
        "heatmap.3(LGT_4heatmap, 
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
}

# Cleanup a few things
`eval \$(/home/ksieber/scripts/change_shell_tmpdir.pl $original_tmpdir)`;
if ( $cluster_fxn eq "pvclust" or $cluster_fxn eq "both" ) { $R->run("stopCluster(cl)"); }
$R->stopR();
print "\n";
print_complete( \%options );

sub read_in_lca_file {
    my $lca_file = shift;

    my $lca_data;    # Data to return

    open( my $IN, "< $lca_file" ) or die "Error: Unable to open the lca_file: $lca_file because: $!\n";
    while (<$IN>) {
        chomp;
        my @split_line = ( split /\t/, $_ );
        my $read_id    = $split_line[0];
        my $lca        = $split_line[$column];
        next if ( !defined $lca or $lca !~ /\w+/ );
        my @split_lca = split( /;/, $lca );
        my @short_lca = ( $lca_cut eq "-1" ) ? @split_lca : splice( @split_lca, 0, $lca_cut );
        my $otu = $short_lca[-1];
        $otu =~ s/\s+/__/g;
        $otu =~ s/\//_/g;
        $otu =~ s/\-/_/g;
        $otu =~ s/\(//g;
        $otu =~ s/\)//g;
        $otu =~ s/\'//g;
        $otu =~ s/\"//g;
        $otu =~ s/\.//g;

        $lca_data->{'otu_by_id'}->{ $links->{$read_id} }->{$otu}++;
        $lca_data->{'total_otu_per_id'}->{ $links->{$read_id} }++;
        $lca_data->{'otus_total'}->{$otu}++;

        $otu_ranks->{$otu} = scalar(@short_lca);    # Global variable. Wrough estimate of specificity of the OTU (higher rank = more specific)
    }
    close $IN;

    return $lca_data;
}

sub merged_micro_lgt_data {
    my $microbiome_data = shift;
    my $lgt_data        = shift;

    my $merged_data;

    # Merge counts of OTUs
    foreach my $otu ( keys %{ $microbiome_data->{'otus_total'} } ) {
        $merged_data->{'otus_total'}->{$otu} = $microbiome_data->{'otus_total'}->{$otu};
    }
    foreach my $otu ( keys %{ $lgt_data->{'otus_total'} } ) {
        $merged_data->{'otus_total'}->{$otu} += $lgt_data->{'otus_total'}->{$otu};
    }

    # Merge ids
    foreach my $id ( keys %{ $microbiome_data->{'otu_by_id'} } ) {
        $merged_data->{'ids'}->{$id} = 1;
    }
    foreach my $id ( keys %{ $lgt_data->{'otu_by_id'} } ) {
        $merged_data->{'ids'}->{$id} = 1;
    }

    return $merged_data;
}

sub reorder_bac_by_counts {
    my $original_top_otu_list = shift;
    my @original_top_otu_list = @{$original_top_otu_list};

    my @new_col_string;
    if ( $sort_bac_by eq "TOTAL" ) {
        @new_col_string = ( sort { ( $merged_data->{otus_total}->{$b} or 0 ) <=> ( $merged_data->{otus_total}->{$a} or 0 ) } @original_top_otu_list );
    }
    elsif ( $sort_bac_by eq "MICRO" ) {
        @new_col_string = ( sort { ( $MICRO_data->{otus_total}->{$b} or 0 ) <=> ( $MICRO_data->{otus_total}->{$a} or 0 ) } @original_top_otu_list );
    }
    elsif ( $sort_bac_by eq "LGT" ) {
        @new_col_string = ( sort { ( $LGT_data->{otus_total}->{$b} or 0 ) <=> ( $LGT_data->{otus_total}->{$a} or 0 ) } @original_top_otu_list );
    }
    return \@new_col_string;
}

sub verify_inputs {
    if ( !defined $options{microbiome_lca} ) { die "Error: Must pass --microbiome_lca=</path/to/microbiome_lca.txt>. Please try again.\n"; }
    if ( !-e $options{microbiome_lca} )      { die "Error: The --microbiome_lca file doesn't exist: $options{microbiome_lca}\n"; }
    if ( !defined $options{lgt_lca} )        { die "Error: Must pass --lgt_lca=</path/to/lgt_lca.txt>. Please try agian.\n"; }
    if ( !-e $options{lgt_lca} )             { die "Error: The --lgt_lca file doesn't exist: $options{lgt_lca}\n"; }
    if ( !defined $options{linker_file} )    { die "Error: Must pass --linker_file=</path/to/file.txt> with Read-ids <tab> participant-id\n"; }
    if ( !-e $options{linker_file} )         { die "Error: The --linker_file doesn't exist: $options{linker_file}\n"; }
    if (    ( defined $options{abundance} )
        and ( $options{abundance} ne 'COUNTS' )
        and ( $options{abundance} ne 'RELATIVE' )
        and ( $options{abundance} ne 'LN' )
        and ( $options{abundance} ne 'LOG' ) )
    {
        die "\n***Error: --abundance has an invalid input: \"$options{abundance}\". Use < COUNTS | RELATIVE | LN | LOG >\n\n";
    }
    if (    ( defined $options{top_otus_by} )
        and ( $options{top_otus_by} ne 'TOTAL' )
        and ( $options{top_otus_by} ne 'BOTH' )
        and ( $options{top_otus_by} ne 'MICRO' )
        and ( $options{top_otus_by} ne 'LGT' ) )
    {
        die "\n***Error: --top_otus_by has an invalid input: \"$options{top_otus_by}\". Use < TOTAL | BOTH | MICRO | LGT >\n\n";
    }
    if (    ( defined $options{sort_bac_by} )
        and ( $options{sort_bac_by} ne 'MICRO' )
        and ( $options{sort_bac_by} ne 'LGT' )
        and ( $options{sort_bac_by} ne 'TOTAL' ) )
    {
        die "\n***Error: --sort_bac_by has an invalid input: \"$options{sort_bac_by}\". Use < TOTAL | MICRO | LGT >\n\n";
    }
    if (    ( defined $options{cluster_fxn} )
        and ( $options{cluster_fxn} ne 'hclust' )
        and ( $options{cluster_fxn} ne 'pvclust' )
        and ( $options{cluster_fxn} ne 'both' ) )
    {
        die "\n***Error: --cluster_fxn has an invalid input: \"$options{cluster_fxn}\". Use < hclust | pvclust | both >\n\n";
    }
    if ( ( defined $options{cluster_by} ) and ( $options{cluster_by} ne 'MICRO' ) and ( $options{cluster_by} ne 'LGT' ) ) {
        die "\n***Error: --cluster_by has an invalid input: \"$options{cluster_by}\". Use < MICRO | LGT >\n\n";
    }
    if ( !defined $options{cluster_num} ) { $options{cluster_num} = "SMART"; }
    if ( ( defined $options{cluster_num} and $options{cluster_num} !~ /^\d+$/ ) and ( $options{cluster_num} ne "SMART" ) ) {
        die "Error: --cluster_num must be numeric or \"SMART\" but was set to: $options{cluster_num}. Please try again.\n";
    }
}

sub initialize_colors_in_R {
    $R->run("reds = rev(rainbow(200))[1:10]");
    $R->run("COLORS = c(rainbow(200)[32:200],reds)");
    $R->run("colorTbl <- rep(NA,0)");
    $R->run("colorTbl[1] <- colors()[90]");
    $R->run("colorTbl[2] <- colors()[553]");
    $R->run("colorTbl[3] <- colors()[47]");
    $R->run("colorTbl[4] <- colors()[12]");
    $R->run("colorTbl[5] <- colors()[461]");
    $R->run("colorTbl[6] <- colors()[526]");
    $R->run("colorTbl[7] <- colors()[429]");
    $R->run("colorTbl[8] <- colors()[234]");
    $R->run("colorTbl[9] <- colors()[652]");
    $R->run("colorTbl[10] <- colors()[624]");
    $R->run("colorTbl[11] <- colors()[550]");
    $R->run("colorTbl[12] <- colors()[173]");
    $R->run("colorTbl[13] <- colors()[103]");
    $R->run("colorTbl[14] <- \"gray\"");
    $R->run("colorTbl[15] <- \"pink\"");
}

sub help {
    print STDERR "\nLGTSeek2heatmaps.PL : This script will take the lca files generated by LGTSeek and create heatmaps for the Microbiome and LGT that have the same id-sorting.\n";
    printf STDERR ( "\t%-30s\t%s\n", "--lgt_lca|l=s",        "/path/to/LGTSeek/lgt_lca.txt" );
    printf STDERR ( "\t%-30s\t%s\n", "--microbiome_lca|m=s", "/path/to/LGTSeek/microbiome_lca.txt" );
    printf STDERR ( "\t%-30s\t%s\n", "--linker_file|L=s",    "/path/to/linker.txt This file has the read-ids and a (participant/analysis/run) id tab-delimted." );
    printf STDERR ( "\t%-30s\t%s\n", "--help_full|?=",       "Full help information." );
    die "\n";
}

sub help_full {
    print STDERR "\nLGTSeek2heatmaps.PL : This script will take the lca files generated by LGTSeek and create heatmaps for the Microbiome and LGT that have the same id-sorting.\n";
    printf STDERR ( "\t%-30s\t%s\n", "--lgt_lca|l=s",        "/path/to/LGTSeek/lgt_lca.txt" );
    printf STDERR ( "\t%-30s\t%s\n", "--microbiome_lca|m=s", "/path/to/LGTSeek/microbiome_lca.txt" );
    printf STDERR ( "\t%-30s\t%s\n", "--linker_file|L=s",    "/path/to/linker.txt - This file has the read-ids and a (participant/analysis/run) ids tab-delimted." );
    printf STDERR (
        "\t%-30s\t%s\n", "--metadata|M=s\@",
        "/path/to/files.txt These metadata files have a line with: \"\#Header <tab> {MetadataName} \" and tab-delimted id : R colors. Accepts multipe -M entries or comma delimted."
    );
    printf STDERR ( "\t%-30s\t%s\n", "--ColorTable=i",     "< 1 | 0 > 1=The Metadata files have numbers to be changed into R colors. 0=Metadata has R color names/codes directly. [1]" );
    printf STDERR ( "\t%-30s\t%s\n", "--abundance|a=s",    "< COUNTS | RELATIVE | LN | LOG > [LN]" );
    printf STDERR ( "\t%-30s\t%s\n", "--top_n_otus=i",     "Top # of otus to keep. [20]" );
    printf STDERR ( "\t%-30s\t%s\n", "--top_otus_by=s",    "< TOTAL | BOTH | MICRO | LGT > Count the top otus based on the  [BOTH]" );
    printf STDERR ( "\t%-30s\t%s\n", "--column|c=i",       "Use lca from this column from the tab-delimted lca-input in BOTH files. [1]" );
    printf STDERR ( "\t%-30s\t%s\n", "--lca_cut=i",        "The lcas are cut short based on this number. [7] " );
    printf STDERR ( "\t%-30s\t%s\n", "--cluster|C=i",      "< 1 | 0 > 1=Cluster the ids based on lca compositions. [0]" );
    printf STDERR ( "\t%-30s\t%s\n", "--cluster_by=s",     "< MICRO | LGT > Cluster the ids by.  [MICRO] " );
    printf STDERR ( "\t%-30s\t%s\n", "--cluster_num|n=i",  "< # > Number of clusters to create OR \"SMART\" (machine learning) [SMART]" );
    printf STDERR ( "\t%-30s\t%s\n", "--cluster_fxn=s",    "< hclust | pvclust > Function to cluster by. [ hclust ] If used, [--cluster=1]." );
    printf STDERR ( "\t%-30s\t%s\n", "--cluster_den=i",    "< 1 | 0 > Display cluster dendrogram. [1]" );
    printf STDERR ( "\t%-30s\t%s\n", "--sort_bac=i",       "< 1 | 0 > 1=Sort the columns by bacterial counts across all data. 0=Semi-Taxonomic rank. [0]" );
    printf STDERR ( "\t%-30s\t%s\n", "--sort_bac_by=s",    "< TOTAL | MICRO | LGT > Data to sort bac by. [TOTAL] If used, [--sort_bac=1]." );
    printf STDERR ( "\t%-30s\t%s\n", "--key=",             "< 1 | 0 > 1= Draw a key for the heatmap. [1] " );
    printf STDERR ( "\t%-30s\t%s\n", "--na_color=s",       "Color to use for lca with NA counts. Takes any \"color\" R:heatmap2 can accept. [#999999]" );
    printf STDERR ( "\t%-30s\t%s\n", "--bottom_margin=s",  "< # > Heatmap bottom_margin.  [16]" );
    printf STDERR ( "\t%-30s\t%s\n", "--right_margin=",    "< # > Heatmap right_margin.   [12]" );
    printf STDERR ( "\t%-30s\t%s\n", "--SidebarColorSize", "< # > Size of the SidebarColors. [2]" );
    printf STDERR ( "\t%-30s\t%s\n", "--output_dir|o=",    "/path/to/directory/for/output [cwd]" );
    printf STDERR ( "\t%-30s\t%s\n", "--output_prefix|p=", "{prefix}_heatmaps.pdf output. [LGTSeek2heatmaps]" );
    printf STDERR ( "\t%-30s\t%s\n", "--help|h",           "Minimal help information and the mandatory inputs." );
    printf STDERR ( "\t%-30s\t%s\n", "--help_full|?=",     "This information." );
    printf STDERR ( "\t%-30s\t%s\n", "--ColorTable_help",  "Help information for the ColorTable." );

    die "\n";
}

sub ColorTable_help {
    die "ColorTable=1 will convert the following numbers from the metadata into the repsective colors.
        1   darkorange
        2   red1
        3   chartreuse
        4   aquamarine4
        5   mediumblue
        6   lightsalmon2
        7   lightseagreen
        8   gray81
        9   yellow
        10  tan4
        11  purple1
        12  gray20
        13  darkseagreen1
        14  gray
        15  pink\n";
}
__END__
