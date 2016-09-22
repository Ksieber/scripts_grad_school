#!/usr/bin/perl -I /home/ksieber/scripts/ -I /home/ksieber/perl5/lib/perl5/

=head1 NAME

bam2regions-of-coverage.pl

=head1 SYNOPSIS

Search an bam for bacterial-human LGT.

=head1 DESCRIPTION


=head1 AUTHORS - Karsten Sieber

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
use POSIX;
use setup_input;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options, 'input|i=s',      'mpileup=s',         'min_cov|c=i', 'make_bams=i', 'sort|S=i', 'threads|t=i', 'sort_mem=s',
    'ref|r=s', 'output_dir|o=s', 'output_prefix|p=s', 'subdirs=i',   'window|w=i',  'help|?'
) or die "Unrecognized command line option. Please try agian.\n";

## Check if the user needs help information
if ( $options{help} ) { &help; }    ## @ end of script
## Check we have the correct input
if ( !$options{input} && !$options{mpileup} ) { die "ERROR: Must pass an input file, use --input=<BAM> or --mpileup=<FILE>\n"; }
if ( $options{input}  && $options{mpileup} )  { die "ERROR: Only pass one input using either --input or --mpileup. Please try again.\n"; }
## Default values
my $lgtseek = LGTSeek->new2( \%options );
my $input = defined $options{input} ? $options{input} : $options{mpileup};
$input =~ s/(\/{2,})/\//g;
my ( $fn, $path, $suffix ) = fileparse( $input, ( @{ $lgtseek->{bam_suffix_list} }, '.mpileup', '.txt' ) );
my $output_prefix = defined $options{output_prefix} ? $options{output_prefix} : "$fn";
my $output_dir    = defined $options{output_dir}    ? $options{output_dir}    : "$path";
my $out           = "$output_dir/$output_prefix";
$out =~ s/(\/{2,})/\//g;
my $p_sort  = defined $options{sort}    ? $options{sort}      : "0";
my $ref     = defined $options{ref}     ? "-f $options{ref} " : undef;
my $min_cov = defined $options{min_cov} ? $options{min_cov}   : "2";
my $window  = defined $options{window}  ? $options{window}    : "0";
my $no_windows_found = 1;

########################
if ( $p_sort == 1 ) {
    run_cmd("samtools sort -@ $lgtseek->{threads} -m $lgtseek->{sort_mem} $input $out\.psort");
    run_cmd("samtools index $out\.psort.bam");
    $input = "$out\.psort.bam";
}
########################
## Open file to find unique regions
my $infh;
if ( $options{mpileup} ) {
    open( $infh, "<", "$input" ) || die "ERROR: Can't open: $options{mpileup}.\n";
}
else {
    open( $infh, "-|", "samtools mpileup $ref -A $input 2>/dev/null" ) || die "ERROR: Can't open: samtools mpileup $ref-A $input\n";
}
########################
# Open output filehandle & print header
my $OFH;
if ( $options{output_dir} || $options{output_prefix} ) {
    open( $OFH, ">", "$out\_regions-of-coverage.txt" ) || die "ERROR: Can't open: $out\_regions-of-coverage.txt.\n";
}
else {
    $OFH = *STDOUT;
}
my $time_stamp = gmtime();
$time_stamp =~ s/\s+/_/g;
print $OFH "## bam2regions-of-coverage.pl MIN_COV=$min_cov WINDOW=$window INPUT=$input TimeStamp=$time_stamp\n";
printf $OFH ( "%-60s%-22s%-24s\n", "## Region", "Mean_Coverage_per_bp", "Total_reads_per_region" );

########################
undef my %start_position;
undef my %previous_position;    ## $previous_position{$chr}=start position of region of coverage
my $final_chr_check;
my @coverage;
my @potential_coverage_extension;

while (<$infh>) {
    chomp;
    my ( $chr, $current_position, $cov ) = ( split /\t/, $_ )[ 0, 1, 3 ];
    $final_chr_check = $chr;
    my $extend_region_position = ( $previous_position{$chr} + 1 + $window );
    if ( $cov >= $min_cov ) {

        # If we move to a diff chr but still cov>min_cov, print and clear to start new window
        if ( !defined $previous_position{$chr} ) {
            undef @potential_coverage_extension;
            &print_region();
            $start_position{$chr}    = $current_position;
            $previous_position{$chr} = $current_position;
            push( @coverage, $cov );
            next;
        }

        # Find the First region and init
        if ( !defined $start_position{$chr} ) {
            push( @coverage, $cov );
            $start_position{$chr}    = $current_position;
            $previous_position{$chr} = $current_position;
            next;
        }

        # If we are moving along the same region (pos+1) add to the region
        if ( $current_position <= $extend_region_position ) {
            push( @coverage, $cov );
            $previous_position{$chr} = $current_position;
            next;
        }

        # If we are in a new region, print the last one and init a new region
        elsif ( $current_position > $extend_region_position ) {
            &print_region;
            $start_position{$chr}    = $current_position;
            $previous_position{$chr} = $current_position;
            push( @coverage, $cov );
        }
    }
    ## If we moved out of a region with coverage
    elsif ( $cov <= $min_cov ) {
        ## If we haven't init. anything yet, keep searching for min_coverage
        next if ( !%start_position or !%previous_position );
        ## If we changed chr, print the region and clear the old window
        if ( !defined $previous_position{$chr} ) {
            undef @potential_coverage_extension;
            &print_region();
            next;
        }
        ## If same chr, but within window size, keep going but don't update the "current position" to maintain window size from last good coverage position
        if ( $current_position < $extend_region_position ) {
            push( @potential_coverage_extension, $cov );
            next;
        }
        else {
            undef @potential_coverage_extension;
            &print_region;
        }
    }
}

#print STDERR "Dumping last region\n";
&print_region;

close $infh;
close $OFH;

if ( $no_windows_found == 1 ) {
    print STDERR "\n\tWarning: No regions with $min_cov where found in: $input\n\n";
}

########################
##### Subroutines ######
########################
sub print_region {
    $no_windows_found = 0;
    push( @coverage, @potential_coverage_extension );
    for my $previous_chr ( keys %previous_position ) {
        my $mean_coverage = floor( Math::NumberCruncher::Mean( \@coverage ) );
        printf $OFH ( "%-60s%-10s", "$previous_chr\:$start_position{$previous_chr}\-$previous_position{$previous_chr}", $mean_coverage );
        if ( defined $options{input} and -e $options{input} ) {
            my $total_read_coverage = run_cmd("samtools view -F0x4 $input \'$previous_chr\:$start_position{$previous_chr}\-$previous_position{$previous_chr}\' | /home/ksieber/bin/WC");
            printf $OFH ( "%-10s", $total_read_coverage );
        }
        print $OFH "\n";
    }
    undef %start_position;
    undef %previous_position;
    undef @coverage;
}

sub help {
    die "    ------------------------------------------------------------
    This script will take a bam file and make a text file with the unique regions of coverage and the average coverage / region.
    Only with input bam, the script will also print total mapped reads to the region.
    OUTPUT goes to STDOUT if no --output_dir or --output_prefix.
    ------------------------------------------------------------
        --input|i=          <BAM>
        --sort=         Sort the input bam.
        --mpileup=          <mpileup of BAM>.
        --ref|r=            Reference for the input.bam. *Recommended*
        --min_cov=          Min. coverage cutoff for regions of interest.
        --window|w=         < # > [0] Number of bases to extend the window. 0=continuous coverge. 
        --make_bams=        Make bams for each region of coverage.
        --output_dir|o=     Directory for output.
         --output_prefix|p=
         --subdirs=
        --help|?
    ------------------------------------------------------------\n";
}

__END__
