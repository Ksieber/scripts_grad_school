#!/usr/bin/perl
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use strict;
use warnings;
use IntervalTree;
use MLDBM 'DB_File';
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;
use mk_dir;
use print_call;
use read_bam;

my %options;
my $results = GetOptions( \%options, 'input|i=s', 'gtf|g=s', 'tree_db|d=s', 'tree_dir=s', 'output|O=s', 'output_dir|o=s', 'output_prefix|p=s', 'help|?' )
    or die "\n*** Error *** Unrecognized command line option. Please try again.\n\n";

if ( $options{help} ) { &help; }
if ( !$options{input} ) { die "Error: Must use --input=< input.bam or file with samtools style regions>. Please try agian or use --help.\n"; }
if ( !$options{gtf}            and !$options{tree_db} )    { die "Error: Must use --gtf=ref.gtf or --tree_db=ref.tree_db Please try again" }
if ( defined $options{gtf}     and !-e $options{gtf} )     { die "Error: This gtf doesn't exist: $options{gtf}\n"; }
if ( defined $options{tree_db} and !-e $options{tree_db} ) { die "Error: This tree_db doesn't exist: $options{tree_db}\n"; }

print_notebook( \%options );

my $input = $options{input};
my ( $in_fn, $in_dir, $in_suf ) = fileparse( $input, qr/\.[^\.]+/ );
my $out_dir = defined $options{output_dir} ? $options{output_dir} : $in_dir;
mk_dir($out_dir);
my $out_pref = defined $options{output_prefix} ? $options{output_prefix} : $in_fn;
my $output   = defined $options{output}        ? $options{output}        : "$out_dir\/$out_pref\_gtf-overlap.txt";
my $OUT_FH;
if ( defined $options{output} or defined $options{output_dir} or defined $options{output_prefix} ) { open( $OUT_FH, ">", "$output" ) or die "Error: Unable to open output: $output\n"; }
else                                                                                               { $OUT_FH = *STDOUT; }

my $tree_db = &load_db();

&analyze_input( $input, $tree_db );

if ( defined $options{output} or defined $options{output_dir} or defined $options{output_prefix} ) { print_complete( \%options ); }

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
########################
##### Subroutines ######
########################

sub load_db {
    my %tree_db;    ## database pointing to each tree

    ## Build the tree_db and trees for the gtf
    if ( defined $options{gtf} and $options{gtf} and !$options{tree_db} ) {
        my $gtf = $options{gtf};
        my ( $gtf_fn, $gtf_dir, $gtf_suff ) = fileparse( $gtf, qr/\.[^\.]+/ );
        $gtf_dir =~ s/\/+$//g;
        my $tree_dir = defined $options{tree_dir} ? $options{tree_dir} : "$gtf_dir\/$gtf_fn\_gtf_overlap_trees";
        $tree_dir =~ s/\/{1,}$//;
        mk_dir($tree_dir);

        # Load from an old ref.tree_db if it exists
        my %tied_db;
        if ( -e "$tree_dir\/$gtf_fn\.tree_db" ) {
            print STDERR "\n*** An existing ref.tree_db ($tree_dir\/$gtf_fn\.tree_db) was found the ref ($gtf) and will be used. ***\n\n";
            tie( %tied_db, "MLDBM", "$tree_dir\/$gtf_fn\.tree_db" );
            %tree_db = %tied_db;
            untie(%tied_db);
            return \%tree_db;
        }

        # Otherwise we will have to create a new ref.tree_db
        else {
            tie( %tied_db, "MLDBM", "$tree_dir\/$gtf_fn\.tree_db" );
        }

        my $tree = IntervalTree->new();
        my $previous_chr;
        my $chr_tree;

        open( IN, "<", $gtf ) or die "Error: Unable to open ref.gtf: $gtf\n";
        while (<IN>) {
            chomp( my $line = $_ );
            next if ( $line =~ /^\#/ );
            my @split_line = split( /\t/, $line );
            my ( $chr, $feature, $start, $end, $attribute ) = ( $split_line[0], $split_line[2], $split_line[3], $split_line[4], $split_line[8] );
            my $hit = (
                {   feature   => $feature,
                    attribute => $attribute
                }
            );

            if ( !$previous_chr ) {
                $previous_chr = $chr;
                $chr_tree     = "$tree_dir\/$gtf_fn\_$chr.tree";
                $tree->addInterval( $hit, $start, $end );
                next;
            }
            elsif ( ( defined $previous_chr ) and ( $chr eq $previous_chr ) ) {
                $tree->addInterval( $hit, $start, $end );
                next;
            }
            elsif ( ( defined $previous_chr ) and ( $chr ne $previous_chr ) ) {

                # Build the previous_chr tree and write to disk
                $tree->buildTree();
                $tree->treeToFile($chr_tree);
                $tree_db{$previous_chr} = $chr_tree;
                undef $tree;

                # Init new tree for the current(next) tree
                $tree = IntervalTree->new();
                $tree->addInterval( $hit, $start, $end );
                $previous_chr = $chr;
                $chr_tree     = "$tree_dir\/$gtf_fn\_$chr.tree";
                next;
            }
            else {
                die "Error: Unsure how to handle this chr: $chr\n";
            }
        }

        # After finish reading the gtf
        # Build the last previous_chr tree and write to disk
        $tree->buildTree();
        $tree->treeToFile($chr_tree);
        $tree_db{$previous_chr} = $chr_tree;

        # Dump the tree_db to disk
        %tied_db = %tree_db;
        untie(%tied_db);

        # Return the tree_db to use as a pointer to the chr_trees
        return \%tree_db;
    }

    elsif ( defined $options{tree_db} and -e $options{tree_db} ) {
        my %tied_db;
        tie( %tied_db, "MLDBM", "$options{tree_db}" );
        %tree_db = %tied_db;
        untie(%tied_db);
        return \%tree_db;
    }
    else {
        die "Error: Unable to determine how to load the ref.gtf tree files.\n";
    }
}

sub analyze_input {
    my $input   = shift;
    my $tree_db = shift;

    my ( $in_fn, $in_dir, $in_suf ) = fileparse( $input, qr/\.[^\.]+/ );
    if ( $in_suf eq ".bam" ) {
        _analyze_bam( $input, $tree_db );
    }
    elsif ( -e $input ) {
        _analyze_txt( $input, $tree_db );
    }
    elsif ( $input =~ /^.+\:\d+\-\d+$/ ) {
        _analyze_region( $input, $tree_db );
    }
    else {
        die "Error: Unable to determine how to analyze the input: $input\n";
    }
}

sub _analyze_txt {
    my $input       = shift;
    my $tree_db_tmp = shift;
    my %tree_db     = %{$tree_db_tmp};
    open( my $in_fh, "<", $input ) or die "Error: Unable to open input.txt: $input\n";
    my $previous_chr;
    my $chr_tree;
    my $tree = IntervalTree->new();
    my $hit;

    while (<$in_fh>) {
        chomp( my $line = $_ );
        next if ( $line =~ /^[\#]+/ );
        my @split_line = split( /\s+/, $line );
        my $id = $split_line[0];
        $id =~ /^([A-Za-z0-9\-\|\_\.]+)\:(\d+)\-(\d+)/;
        my ( $chr, $start, $stop ) = ( $1, $2, $3 );
        if ( $chr =~ /gi\|\d+\|\w{2,3}\|([A-Za-z0-9]+)\.?\d+\|/ ) { $chr = $1; }

        if ( !$previous_chr ) {
            if ( !$tree_db{$chr} ) {
                printf $OUT_FH ( "%-60s%-20s%-40s", $id, "NaN", "No $chr in ref.gtf" );
                print $OUT_FH "\n";
                next;
            }
            $previous_chr = $chr;
            $chr_tree     = $tree_db{$chr};
            $tree->treeFromFile($chr_tree);
            my @hits = $tree->searchInterval( $start, $stop );
            if (@hits) {
                foreach my $hit (@hits) {
                    printf $OUT_FH ( "%-60s%-20s\t%s", $id, $hit->[2]->{feature}, $hit->[2]->{attribute} );    ## %-60s%-40s%-300s
                    print $OUT_FH "\n";
                    next;
                }
            }
        }
        elsif ( ( defined $previous_chr ) and ( $chr eq $previous_chr ) ) {
            my @hits = $tree->searchInterval( $start, $stop );
            if (@hits) {
                foreach my $hit (@hits) {
                    printf $OUT_FH ( "%-60s%-20s\t%s", $id, $hit->[2]->{feature}, $hit->[2]->{attribute} );
                    print $OUT_FH "\n";
                    next;
                }
            }
        }
        elsif ( ( defined $previous_chr ) and ( $chr ne $previous_chr ) ) {
            if ( !$tree_db{$chr} ) {
                printf $OUT_FH ( "%-60s%-20s%-40s", $id, "NaN", "No $chr in ref.gtf" );
                print $OUT_FH "\n";
                next;
            }
            $previous_chr = $chr;
            $tree         = IntervalTree->new();
            $chr_tree     = $tree_db{$chr};
            $tree->treeFromFile($chr_tree);
            my @hits = $tree->searchInterval( $start, $stop );
            if (@hits) {
                foreach my $hit (@hits) {
                    printf $OUT_FH ( "%-60s%-20s\t%s", $id, $hit->[2]->{feature}, $hit->[2]->{attribute} );
                    print $OUT_FH "\n";
                    next;
                }
            }
        }
        else {
            die "Error: Unable to determine how to do things.\n";
        }
    }
    close $in_fh;
}

sub _analyze_bam {
    my $bam         = shift;
    my $tree_db_tmp = shift;
    my %tree_db     = %{$tree_db_tmp};
    my $bam_fh      = open_bam($bam);
    my $previous_chr;
    my $chr_tree;
    my $tree = IntervalTree->new();

    while ( my $read = read_bam($bam_fh) ) {
        next if ( $read->{flag}->{qunmapped} );
        my $stop = $read->{position} + length( $read->{sequence} );
        if ( !$previous_chr ) {
            if ( !$tree_db{ $read->{chr} } ) {
                printf $OUT_FH ( "%-60s%-20s\t%s", $read->{id}, "NaN", "No $read->{chr} in ref.gtf" );
                print $OUT_FH "\n";
                next;
            }
            $previous_chr = $read->{chr};
            $chr_tree     = $tree_db{ $read->{chr} };
            $tree->treeFromFile($chr_tree);
            my @hits = $tree->searchInterval( $read->{position}, $stop );
            if (@hits) {
                foreach my $hit (@hits) {
                    printf $OUT_FH ( "%-60s%-36s%-20s\t%s", "$read->{chr}\:$read->{position}\-$stop", $read->{id}, $hit->[2]->{feature}, $hit->[2]->{attribute} );
                    print $OUT_FH "\n";
                    next;
                }
            }
        }
        elsif ( ( defined $previous_chr ) and ( $read->{chr} eq $previous_chr ) ) {
            my @hits = $tree->searchInterval( $read->{position}, $stop );
            if (@hits) {
                foreach my $hit (@hits) {
                    printf $OUT_FH ( "%-60s%-36s%-20s\t%s", "$read->{chr}\:$read->{position}\-$stop", $read->{id}, $hit->[2]->{feature}, $hit->[2]->{attribute} );
                    print $OUT_FH "\n";
                    next;
                }
            }
        }
        elsif ( ( defined $previous_chr ) and ( $read->{chr} ne $previous_chr ) ) {
            if ( !$tree_db{ $read->{chr} } ) {
                printf $OUT_FH ( "%-60s%-36s%-20s", $read->{id}, "NaN", "No $read->{chr} in ref.gtf" );
                print $OUT_FH "\n";
                undef $previous_chr;
                next;
            }
            $previous_chr = $read->{chr};
            $chr_tree     = $tree_db{ $read->{chr} };
            $tree         = IntervalTree->new();
            $tree->treeFromFile($chr_tree);
            my @hits = $tree->searchInterval( $read->{position}, $stop );
            if (@hits) {
                foreach my $hit (@hits) {
                    printf $OUT_FH ( "%-60s%-36s%-20s\t%s", "$read->{chr}\:$read->{position}\-$stop", $read->{id}, $hit->[2]->{feature}, $hit->[2]->{attribute} );
                    print $OUT_FH "\n";
                    next;
                }
            }
        }
        else {
            die "Error: Unable to determine how to do things.\n";
        }
    }
    close $bam_fh;
}

sub _analyze_region {
    my $region      = shift;
    my $tree_db_tmp = shift;
    my %tree_db     = %{$tree_db_tmp};

    $region =~ /^([A-Za-z0-9\-\|\_\.]+)\:(\d+)\-(\d+)/;
    my ( $chr, $start, $stop ) = ( $1, $2, $3 );
    if ( $chr =~ /gi\|\d+\|\w{2,3}\|([A-Za-z0-9]+)\.?\d+\|/ ) { $chr = $1; }

    my $tree     = IntervalTree->new();
    my $chr_tree = $tree_db{$chr};
    $tree->treeFromFile($chr_tree);

    my @hits = $tree->searchInterval( $start, $stop );
    if (@hits) {
        foreach my $hit (@hits) {
            printf $OUT_FH ( "%-60s%-20s\t%s", $region, $hit->[2]->{feature}, $hit->[2]->{attribute} );
            print $OUT_FH "\n";
            next;
        }
    }
    else {
        print $OUT_FH "%-60s%-20s\t%s", $region, "NaN", "No hits for $region in ref.gtf\n";
    }
}

sub help {
    die "\nThis script will take an input file and find the overlap regions from a ref.gtf annotation file. 
    --input=            Input can either be bam or file with samtools style regions in the first coloumn, 1 per line. If the input is a bam, it *must* be position sorted.
    --gtf|g=            A ref annotated.gtf file.
    --tree_db=          A database (ref.tree_db) will be created for a ref.gtf file the first time running the script. After that, the script can be pointed to/will find the ref.tree_db and run faster.
    --tree_dir=         Directory to create the ref.tree_db in. Default = ref_dir
    --output_dir=       Directory for output.
    --output_prefix=    output_prefix.txt\n\n"
}
