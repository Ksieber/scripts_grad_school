#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/

=head1 NAME

blast2lca.pl

=head1 SYNOPSIS

Calculate the LCA from hits in a blast-m8 report.

=head1 DESCRIPTION


=head1 AUTHORS - Karsten Sieber

e-mail: Karsten.sieber@gmail.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with "_"

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,
    'input|i=s',    # Comma separated list of files
    'input_list|I=s',
    'evalue_cutoff|e=s',
    'best_hits_only=i',
    'combine_PE_lca=i',
    'Qsub|q=i',
    'project=s',
    'sub_name=s',
    'sub_mail=s',
    'print_hostname|ph=i',
    'output_dir|o=s',
    'output_prefix=s',
    'subdirs=i',
    'taxon_host=s',
    'taxon_dir=s',
    'taxon_idx_dir=s',
    'conf_file=s',
    'verbose=i',
    'help|?',
) or die "Error: Unrecognized command line option. Please try again.\n";

use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/', '/local/projects-t3/HLGT/scripts/lgtseek/lib/', '/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/' );
use print_call;
print_hostname( \%options );    ## This is useful for trouble shooting grid nodes that might be missing modules for LGTSeek etc.

use LGTSeek;
use run_cmd;
use setup_input;
use Carp;
$Carp::MaxArgLen = 0;
use File::Basename;

if ( $options{help} ) { &help; }    ## @ end of script
if ( !$options{input} && !$options{input_list} ) { die "Error: Please give an input blast-m8 file with --input=<FILE> or --input_list=<LIST>. Try again or use --help_full.\n"; }

## Qsub:
## Initialize LGTSeek.pm
my $lgtseek = LGTSeek->new2( \%options );
if ( $options{Qsub} ) {
    if ( !$options{sub_name} ) { $options{sub_name} = "m8_2lca"; }
    if ( !$options{project} )  { $options{project}  = $lgtseek->{project}; }
    Qsub_script( \%options );
}

## Print the script call
print_call( \%options );

## Setup array ref of inputs
my $inputs = setup_input( \%options );

# Setup a few default values
my $evalue_cutoff  = $options{evalue_cutoff}  ? $options{evalue_cutoff}  : "1";
my $bho            = $options{best_hits_only} ? $options{best_hits_only} : "0";
my $combine_PE_lca = $options{combine_PE_lca} ? $options{combine_PE_lca} : "0";

# Calculate LCA's
foreach my $input (@$inputs) {
    my ( $name, $path, $suf ) = fileparse( $input, $lgtseek->{suffix_regex} );
    chomp $name;

    ## Setup output directory
    if ( !$lgtseek->{output_dir} ) { $lgtseek->{output_dir} = "$path/lgtseq/"; }
    if ( $lgtseek->{subdirs} ) {
        $lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}");
        $lgtseek->{output_dir} = "$options{output_dir}/" . "$name/";
        $lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}");
    }

    my $lca = $lgtseek->blast2lca(
        {   blast          => $input,
            output_dir     => $lgtseek->{output_dir},
            output_prefix  => $options{output_prefix},
            evalue_cutoff  => $evalue_cutoff,
            best_hits_only => $bho,
            combine_PE_lca => $combine_PE_lca,
            verbose        => $lgtseek->{verbose},
        }
    );
}

print_complete( \%options );

sub help {
    die "Help: This script will take a blast-m8 report and calculate the LCA for each unique ID.
    ----------------------------------------------------------------------------------------------------
        --input|i=              <Input> Accepts blast -m8 reports. 
        --input_list=           List of blast -m8 reports, 1 file per line.
    ----------------------------------------------------------------------------------------------------
        --evalue_cutoff|e=      # For the highest evalue allowed.
        --best_hits_only=       <0|1> [0] 1= Parse the blast hit for only best hits.
        --combine_PE_lca=       <0|1> [0] 1= Merge both PE reads into 1 LCA with the conservative and liberal methods.
    ----------------------------------------------------------------------------------------------------
        --output_dir|o=         Directory for all output. Will be created if it doesn't exist.
          --output_prefix=      Name of prefix for the output. 
          --subdirs=            <0|1> [0] 1 = Create a new subdirectory for each input in the output directory.
    ----------------------------------------------------------------------------------------------------
        --Qsub|q=               <0|1> [0] 1 = Submit the job to the grid.
          --sub_name=           Name of the job to submit.
          --project=            [jdhotopp-lab] Grid project to use.
          --sub_mail=           [0] 1= email user\@som.umaryland.edu when job is complete & with stats. Can also specify --sub_mail=specific\@email.foo
    ----------------------------------------------------------------------------------------------------    
        --conf_file=            [/home/USER/.lgtseek.conf] 
        --verbose|V             <0|1> 1=Verbose reporting of progress. 0=Turns off most reports.
        --help|?
        Note: Make sure to have ~/.lgtseek.conf taxon files the same as the blast db used.
    ----------------------------------------------------------------------------------------------------\n";
}

