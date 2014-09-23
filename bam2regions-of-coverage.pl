#!/usr/bin/perl -I /home/ksieber/scripts/ -I /home/ksieber/perl5/lib/perl5/

=head1 NAME

bam2regions-of-coverage.pl

=head1 SYNOPSIS

Search an bam for bacterial-human LGT.

=head1 DESCRIPTION


=head1 AUTHORS - Karsten Sieber & David Riley

e-mail: ksieber@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with "_"

=cut

use warnings;
no warnings 'uninitialized';
use strict;
use Math::NumberCruncher;
use lib qw(/local/projects-t3/HLGT/scripts/lgtseek/lib/ /opt/lgtseek/lib/);    ### May need to change this depending on where the script is being run
use LGTSeek;
use Time::SoFar;
use run_cmd;
use setup_input;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'mpileup=s', 'min_cov|c=i', 'make_bams=i', 'pos_sort|S=i', 'threads|t=i', 'sort_mem=s', 'ref|r=s', 'output_dir|o=s', 'output_prefix|p=s', 'subdirs=i', 'help|?' )
    or die "Unrecognized command line option. Please try agian.\n";

## Check if the user needs help information
if ( $options{help} ) { &help; }                                               ## @ end of script
## Check we have the correct input
if ( !$options{input} && !$options{mpileup} ) { die "ERROR: Must pass an input file, use --input=<BAM> or --mpileup=<FILE>\n"; }
if ( $options{input}  && $options{mpileup} )  { die "ERROR: Only pass one input using either --input or --mpileup. Please try again.\n"; }
## Default values
my $lgtseek = LGTSeek->new2( \%options );
my $input = defined $options{input} ? $options{input} : $options{mpileup};
my ( $fn, $path, $suffix ) = fileparse( $input, ( @{ $lgtseek->{bam_suffix_list} }, '.mpileup', '.txt' ) );
my $output_prefix = defined $options{output_prefix} ? $options{output_prefix} : "$fn";
my $output_dir    = defined $options{output_dir}    ? $options{output_dir}    : "$path";
my $out           = "$output_dir/$output_prefix";
my $p_sort        = defined $options{pos_sort}      ? $options{pos_sort}      : "0";
my $ref           = defined $options{ref}           ? "-f $options{ref} "     : undef;
my $min_cov       = defined $options{min_cov}       ? $options{min_cov}       : "2";

########################
if ( $p_sort == 1 ) {
    run_cmd("samtools sort -@ $lgtseek->{threads} -m $lgtseek->{sort_mem} $input $out\.pos-sort");
    run_cmd("samtools index $out\.pos-sort.bam $out\.pos-sort.bai");
    $input = "$out\.pos-sort.bam";
}
########################
## Open file to find unique regions
my $infh;
if ( $options{mpileup} ) {
    open( $infh, "<", "$input" ) || die "ERROR: Can't open: $options{mpileup}.\n";
}
else {
    open( $infh, "-|", "samtools mpileup $ref -A $input" ) || die "ERROR: Can't open: samtools mpileup $ref-A $input\n";
}
######################## 
my $OFH;
if ( $options{output_dir} || $options{output_prefix} ) {
    open( $OFH, ">", "$out\_regions-of-coverage.txt" ) || die "ERROR: Can't open: $out\_regions-of-coverage.txt.\n";
}
else {
    $OFH = *STDOUT;
}
########################
undef my %start_position;
undef my %previous_position;    ## $previous_position{$chr}=start position of region of coverage
my $final_chr_check;
my @coverage;
while (<$infh>) {
    chomp;
    my ( $chr, $current_position, $cov ) = ( split /\t/, $_ )[ 0, 1, 3 ];
    $final_chr_check = $chr;
    if ( $cov >= $min_cov ) {
        ## If we move to a diff chr but still cov>min_cov, print and clear to start new
        my $new_chr = 0;
        for my $previous_chr ( keys %previous_position ) {
            if ( $chr !~ /$previous_chr/ ) {

                #print STDERR "New Chr: $_\n";
                &print_region();
                $start_position{$chr}    = $current_position;
                $previous_position{$chr} = $current_position;
                push( @coverage, $cov );
                $new_chr = 1;
            }
        }
        next if ( $new_chr == 1 );
        ## Find the First region and init
        if ( !$start_position{$chr} ) {

            #print STDERR "Init start of a region: $_\n";
            push( @coverage, $cov );
            $start_position{$chr}    = $current_position;
            $previous_position{$chr} = $current_position;
            next;
        }
        ## If we are moving along the same region (pos+1) add to the region
        my $extend_region_position = ( $previous_position{$chr} + 1 );
        if ( $current_position == $extend_region_position ) {

            #print STDERR "Adding to region: $_\n";
            push( @coverage, $cov );
            $previous_position{$chr} = $current_position;
            next;
        }
        ## If we are in a new region, print the last one and init a new region
        if ( $current_position != $extend_region_position ) {

            #print STDERR "New region: $_\n";
            &print_region;
            $start_position{$chr}    = $current_position;
            $previous_position{$chr} = $current_position;
            push( @coverage, $cov );
        }
        ## If we moved out of a cov. region, print it and clear regions
    }
    elsif ( $cov <= $min_cov ) {
        next if ( !%start_position || !%previous_position );

        #print STDERR "End of region: $_\n";
        &print_region;
    }
}

#print STDERR "Dumping last region\n";
&print_region;

close $infh;
close $OFH;
########################
##### Subroutines ######
########################
sub print_region {
    for my $previous_chr ( keys %previous_position ) {
        my $mean_coverage = Math::NumberCruncher::Mean( \@coverage );
        print $OFH "$previous_chr\:$start_position{$previous_chr}\-$previous_position{$previous_chr}\t$mean_coverage\n";
    }
    undef %start_position;
    undef %previous_position;
    undef @coverage;
}

sub help {
    die "    ------------------------------------------------------------
    This script will take a bam file and make a text file with the unique regions of coverage.
    OUTPUT goes to STDOUT if no --output_dir or --output_prefix.
    ------------------------------------------------------------
        --input=            <BAM>
        --pos_sort=         Sort the input bam.
        --mpileup=          <mpileup of BAM>.
        --min_cov=          Min. coverage cutoff for regions of interest.
        --make_bams=        Make bams for each region of coverage.
        --output_dir=       Directory for output.
         --output_prefix=
         --subdirs=
        --help
    ------------------------------------------------------------\n";
}

__END__
