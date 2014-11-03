#!/usr/local/bin/perl -w
use strict;
use run_cmd;
use mk_dir;
use setup_input;
use print_call;
use List::MoreUtils qw( minmax );
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'input_list|I=s', 'output_dir|o=s', 'output_prefix|p=s', 'min_max=s', 'chr_min_max=s', 'help|?', )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) {
    die "HELP: This script will take a list of mpileup files and merge them by position.
	--input|i=		    Directory to find all *.mpileup files from.
	--input_list|I=		List of the files to merge. 1 file per line. \tREQUIRED.
	--output_dir|o=		Directory to put the output.\t\t\tOPTIONAL.
	--output_prefix|p=	Prefix for output. \$prefix\_merged-mpileup.txt \t\t\t\tOPTIONAL.
	--min_max=		<100-200> Override the min/max values for ALL chr (warning, be smart with this).
	--chr_min_max=		<chr:100-200> Override the min/max values for ONE chr.
    --help|?\n";
}

if ( !$options{input_list} && !$options{input} ) {
    die "Must give an input with --input=<DIR> --input_list=<LIST of mpileups>\n";
}

my @input_list;
if ( $options{input} ) {
    chomp( @input_list = `find $options{input} -name \'*.mpileup\'` );
}
elsif ( $options{input_list} ) {
    open( LIST, "<", $options{input_list} ) or die "Unable to open the input list: $options{input_list} because: $!\n";
    while (<LIST>) {
        chomp;
        push( @input_list, $_ );
    }
}
my ( $foo, $bar, $salad ) = fileparse( $input_list[0], 'mpileup' );
my $output_dir = $options{output_dir} ? $options{output_dir} : "$bar";
mk_dir($output_dir);
my $out_pref = $options{output_prefix} ? $options{output_prefix} : "$foo";
my $output = "$output_dir\/$out_pref\_merged-mpileup.txt";

###############################
## May need to adjust the split for "different" mpileup's.
my %foo;
my %chr_min;
my %chr_max;
my @header_list;

foreach my $input (@input_list) {
    my $header;
    if ( $input =~ /(hg19_SRR\d+)/ ) {
        $header = $1;
    }
    else {
        my @suffix_list = ( '.mpileup', '.COVERAGE.txt', 'FF_mpileup', 'RR_mpileup' );    ## May need to change this depending on input
        my ( $file, $path, $suffix ) = fileparse( $input, @suffix_list );
        $header = $file;
    }
    push( @header_list, $header );
    open( my $handle, "<", "$input" ) or die "Unable to open: $input because: $!\n";
    while (<$handle>) {
        chomp;
        if ( $_ !~ /\w*/ ) {
            $foo{0}->{$header} = 0;
        }
        my ( $chromosome, $bp, $coverage ) = ( split( /\t/, $_ ) )[ 0, 1, 3 ];
        my $position = "$chromosome;$bp";                                                 ## pos = chr ; bp
        $foo{$position}->{$header} = $coverage;  ## May need to change this line if the mpileup is a hack. =$f[coverage col.] $f[3] = default mpileup [Hash{chr:bp} -> sample_header = mpileup coverage]
        if ( !$chr_min{$chromosome} ) {
            $chr_min{$chromosome} = $bp;
        }
        if ( !$chr_max{$chromosome} ) {
            $chr_max{$chromosome} = $bp;
        }
        if ( $bp < $chr_min{$chromosome} ) {
            $chr_min{$chromosome} = $bp;
        }
        if ( $bp > $chr_max{$chromosome} ) {
            $chr_max{$chromosome} = $bp;
        }
    }
}

my $override_chr;
my $override_min;
my $override_max;

if ( $options{min_max} ) {
    $options{min_max} =~ /(\d+)\-(\d+)$/;
    $override_min = $1;
    $override_max = $2;
    foreach my $keys ( keys %chr_min ) {
        $chr_min{$keys} = $override_min;
    }
    foreach my $keys ( keys %chr_max ) {
        $chr_max{$keys} = $override_max;
    }
}

if ( $options{chr_min_max} ) {
    $options{chr_min_max} =~ /^([A-Za-z0-9_\.\|]+)\:(\d+)\-(\d+)$/;
    $override_chr = $1;
    $override_min = $2;
    $override_max = $3;
    if ( !$chr_min{$override_chr} ) { die "ERROR: The chr you specified wasn't seen.\n"; }
    $chr_min{$override_chr} = $override_min;
    $chr_max{$override_chr} = $override_max;

    # print STDERR "$override_chr\t$chr_min{$override_chr}\t$chr_max{$override_chr}\n";

}

#foreach my $k (keys %chr_min){
#    print STDERR "MIN:$chr_min{$k}\tMAX:$chr_max{$k}\n";
#}

open( OUT, ">", "$output" ) or die "Unable to open output: $output because: $1\n";
print OUT "Chr;Position";
foreach my $header2 (@header_list) {
    print OUT "\t$header2";
}
print OUT "\n";

foreach my $chr ( sort keys %chr_min ) {
    for ( my $i = $chr_min{$chr}; $i < $chr_max{$chr} + 1; $i++ ) {
        my $pos = "$chr;$i";
        print OUT "$pos\t";
        foreach my $header3 (@header_list) {
            if ( $foo{$pos}->{$header3} ) {
                print OUT "\t$foo{$pos}->{$header3}";
            }
            else {
                print OUT "\t0";
            }
        }
        print OUT "\n";
    }
}
print_complete( \%options );
