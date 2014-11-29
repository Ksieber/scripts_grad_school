#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/

=head1 NAME

bam_prep_for_lgt_seq.pl

=head1 SYNOPSIS

Split, encrypt, and/or prelim filter a bam for lgt_seq.pl 

=head1 DESCRIPTION

Split, encrypt, and/or prelim filter a bam for lgt_seq.pl 

=head1 AUTHOR - Karsten Sieber

e-mail: ksieber@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

use warnings;
no warnings 'uninitialized';
use strict;
use Carp;
$Carp::MaxArgLen = 0;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use print_call;
use lib ( "/local/projects-t3/HLGT/scripts/lgtseek/lib/", "/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/" );
use LGTSeek;
use setup_input;
use run_cmd;
use print_call;

my %options;
my $results = GetOptions( \%options, 'input=s', 'input_list=s', 'output_dir=s', 'threads=i', 'rate_limit=s', 'cghub_key=s', 'Qsub=i', 'sub_mail=s', 'help|?', 'verbose=i', 'Qsub_iterate=i', )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) {
    die "This script will download the data from a cgquery.xml file.
    --input=                Either an analysis_id or xml file from cgquery.
    --input_list=           List of analysis_id's, 1 id / line.
    --output_dir=           /dir/for/output/
    --threads=              Number of threads to use for download
    --Qsub=                 1= Submit job to the grid.
      --sub_name=
    --rate_limit        
    --cghub_key=            [~/.lgtseq.config]
    --verbose=              <0|1> [0] 1= Verbose reporting.
    --help\n";
}

if ( !$options{input} && !$options{input_list} ) { die "Must use --input\n"; }
if ( !$options{output_dir} ) { die "Must give use --output_dir=\"/dir/for/output/\"\n"; }
if ( $options{Qsub} or $options{Qsub_iterate} ) { Qsub_script( \%options ) }

my $lgtseq = LGTSeek->new2( \%options );

print_call( \%options );

if ( -e $options{input} && $options{input} =~ /\.xml$/ ) {
    my $bams = $lgtseq->downloadCGHub(
        {   xml        => $options{input},
            output_dir => $options{output_dir},
            threads    => $lgtseq->{threads},
            rate_limit => $lgtseq->{rate_limit},
            cghub_key  => $lgtseq->{cghub_key},
        }
    );
}
elsif ( $options{input} ) {
    my $bams = $lgtseq->downloadCGHub(
        {   analysis_id => $options{input},
            output_dir  => $options{output_dir},
            threads     => $lgtseq->{threads},
            rate_limit  => $lgtseq->{rate_limit},
            cghub_key   => $lgtseq->{cghub_key},
        }
    );
}
elsif ( $options{input_list} && -e $options{input_list} ) {
    my $input = setup_input( \%options );
    foreach my $analysis_id (@$input) {
        my $bams = $lgtseq->downloadCGHub(
            {   analysis_id => $analysis_id,
                output_dir  => $options{output_dir},
                threads     => $lgtseq->{threads},
                rate_limit  => $lgtseq->{rate_limit},
                cghub_key   => $lgtseq->{cghub_key},
            }
        );
    }
}
else {
    confess "Error: Unable to properly process input.\n";
}

print_complete( \%options );

