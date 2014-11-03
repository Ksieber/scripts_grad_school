#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/

=head1 NAME

lgt_hybrid_finder.pl

=head1 SYNOPSIS

Search for reads ON the integration point. ie Half Donor Half Host.

=head1 DESCRIPTION


=head1 AUTHOR - Karsten B. Sieber

e-mail: ksieber@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,
    'input|i=s',    ## Comma separated list of files
    'output_dir|o=s',
    'prinseq_filter=i',
    'bin_dir=s',
    'samtools_bin=s',
    'ergatis_bin=s',
    'prinseq_bin=s',
    'donor_lineage=s',
    'host_lineage=s',
    'taxon_host=s',
    'taxon_dir=s',
    'taxon_idx_dir=s',
    'path_to_blastdb=s',
    'Qsub=i',
    'sub_mail=s',
    'conf_file=s',
    'print_hostname|ph=i',
    'help|h'
) or die "Error: Unrecognized command line option. Please try again.\n";

use print_call;
$options{print_hostname} = $options{print_hostname} ? $options{print_hostname} : 0;
print_hostname( \%options );    ## This is useful for trouble shooting grid nodes that might be missing modules for LGTSeek etc.

use lib ( "/local/projects-t3/HLGT/scripts/lgtseek/lib/", "/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/" );
use run_cmd;
use LGTSeek;
use File::Basename;

if ( $options{help} ) { &help(); }

#if(!$options{input}){die "Error: Please give an input with --input=<.fa or .bam>. Try again or use --help.\n";} ## KBS 01.17.14

if ( $options{Qsub} == 1 ) { Qsub_script( \%options ) }
print_call( \%options );

# Create an lgtseek object
my $lgtseek = LGTSeek->new2( \%options );

# Take care of the inputs
## Setup Default paths for references and bins:
my $prinseq_filter = defined $options{prinseq_filter} ? $options{prinseq_filter} : "0";
my ( $name, $path, $suf ) = fileparse( $options{input}, $lgtseek->{suffix_regex} );

my $input;
print STDERR "=====  BESTBLAST2  =====\n";
## If we have a bam input, make a fasta for blast filterting
if ( $suf eq '.bam' ) {
    if ($prinseq_filter) {
        my $filtered_bam = $lgtseek->prinseqFilterBam(
            {   input_bam  => $options{input},
                output_dir => "$lgtseek->{output_dir}/prinseq_filter/",
            }
        );
        $input = $lgtseek->sam2Fasta(
            {   input      => "$filtered_bam->{file}",
                output_dir => "$lgtseek->{output_dir}/best_blast/",
            }
        );
    }
    else {
        $input = $lgtseek->sam2Fasta(
            {   input      => $options{input},
                output_dir => "$lgtseek->{output_dir}/best_blast/",
            }
        );
    }
}    ## If we have a fasta input, set for blast filtering
elsif ( $suf eq '.fa' || '.fasta' ) {
    $input = $options{input};
}    ## If we have a blast_list_file, dont make a fasta file
elsif ( $suf ne '.list' ) {
    die "Error: Could not determine input file type. Please use a: .bam .fa or blast_list_file\n";
}

my $best_blasts;
if ( $suf ne '.list' ) {

    # Blast & get best hit
    $best_blasts = $lgtseek->bestBlast2(    ## KBS 01.17.14
        {   db         => $lgtseek->{path_to_blastdb},
            lineage1   => $lgtseek->{donor_lineage},             ## 'Bacteria',
            lineage2   => $lgtseek->{host_lineage},              ## 'Eukaryota',
            fasta      => $input,
            output_dir => "$lgtseek->{output_dir}/best_blast/"
        }
    );
}
elsif ( $suf eq '.list' ) {
    $best_blasts->{list_file} = $options{input};
}

# Now run lgtfinder
print STDERR "===== LGTFINDER =====\n";
my $lgt_blast_validated = $lgtseek->runLgtFinder(
    {   lineage1        => $lgtseek->{donor_lineage},              ## 'Eukaryota',
        lineage2        => $lgtseek->{host_lineage},               ## 'Bacteria',
        input_file_list => $best_blasts->{list_file},
        output_prefix   => $name,
        max_overlap     => $lgtseek->{max_overlap},
        min_length      => $lgtseek->{min_length},
        output_dir      => "$lgtseek->{output_dir}/lgt_finder/",
    }
);

if ( $suf eq '.bam' ) {
    my $validated_lgts = $lgtseek->validated_bam(
        {   input    => $lgtseek->{input},
            by_clone => $lgt_blast_validated->{by_clone},
            output   => "$lgtseek->{output_dir}/$name\_lgt_valid_blast.bam"
        }
    );
}

print_complete( \%options );

sub help {
    die "Help: This script will takes a bam and identifies bacterial human LGT.
        --input|i=              < BAM, FASTA, or blast_list_file >
        --output_dir|o=             </dir/for/output/> [cwd]
        --prinseq_filter=           <0|1> [0] 1=prinseq_filter input bam (not implemented for fasta yet)
        --Qsub=                     <0|1> [0] 1=Submit job to SGE grid.
         --sub_mail=             [0] 1= email user\@som.umaryland.edu when job is complete & with stats. Can also specify --sub_mail=specific\@email.foo 
        --taxon_host=
        --taxon_dir=
        --taxon_idx_dir=
        --path_to_blastdb=
        --bin_dir=
        --samtools_bin=
        --ergatis_bin=
        --conf_file=                [~/.lgtseek.conf]
        --help\n";
}

__END__
