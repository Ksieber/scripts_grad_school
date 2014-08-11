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
    \%options,             'fasta|f=s',            'bam|b=s',        'db=s',               'formatdb=i', 'm8=i',
    'evalue_cutoff=s',     'best_hits_only|bho=i', 'output_dir|o=s', 'output_prefix|op=s', 'Qsub|Q=i',   'sub_name=s',
    'print_hostname|ph=i', 'subdirs=i',            'threads|t=i',    'conf_file=s',        'verbose=i',  'help|h',
) or die "Error: Unrecognized command line option. Please try again.\n";

use lib ( "/local/projects-t3/HLGT/scripts/lgtseek/lib/", "/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/" );
use print_call;
if ( $options{ph} ) { print_hostname( \%options ); }    ## This is useful for trouble shooting grid nodes that might be missing modules for LGTSeek etc.

use LGTSeek;
use run_cmd;
use setup_input;
use Carp;
$Carp::MaxArgLen = 0;
use File::Basename;
use Bio::SearchIO;
use Bio::SearchIO::blasttable;
use Bio::SearchIO::Writer::ResultTableWriter;

if ( $options{help} ) { &help; }                        ## @ end of script
if ( !$options{fasta} && !$options{bam} ) { confess "Error: Please give an input with --fasta=<FASTA> or --bam=<BAM>.\n"; }

## Print the script call
print_call( \%options );
if ( $options{Qsub} ) {
    if ( !$options{sub_name} ) { $options{sub_name} = "blastn"; }
    Qsub_script( \%options );
}

## Initialize LGTSeek.pm
my $lgtseek = LGTSeek->new2( \%options );
$lgtseek->{verbose} = 1;

## Setup output directory
my $input = defined $options{bam} ? $options{bam} : $options{fasta};
my ( $name, $path, $suf ) = fileparse( $input, $lgtseek->{suffix_regex} );
chomp $name;
if ( !$lgtseek->{output_dir} ) { $lgtseek->{output_dir} = "$path/blastn/"; }
$lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}");
if ( $lgtseek->{subdirs} ) {
    $lgtseek->{output_dir} = "$options{output_dir}/" . "$name/";
    $lgtseek->_run_cmd("mkdir -p $lgtseek->{output_dir}");
}
print_notebook( \%options );

my $output_prefix = defined $options{output_prefix} ? $options{output_prefix} : "$name\_blastn";
if ( $options{formatdb} ) {
    run_cmd("formatdb -p F -i $options{db}");
}

# Convert Bam -> Fasta for Blast
my $fasta = defined $options{bam} ? $lgtseek->sam2fasta( { input => $options{bam}, output_dir => $lgtseek->{output_dir} } ) : $options{fasta};

my $evalue_cutoff = defined $options{evalue_cutoff}  ? $options{evalue_cutoff} : '1';
my $m8            = defined $options{m8}             ? $options{m8}            : '0';
my $bho           = defined $options{best_hits_only} ? '1'                     : '0';

if ( $m8 == 1 ) {

    # Run Blast
    open( IN, "blastall -p blastn -a $lgtseek->{threads} -e $evalue_cutoff -T F -d $options{db} -m8 -i $fasta |" )
        || confess "Error: Unable to start blast fh: blastall -p blastn -a $lgtseek->{threads} -e 1 -T F -d $options{db} -m8 -i $fasta\n";
    open( OUT, ">$lgtseek->{output_dir}/$output_prefix\-m8.txt" ) or confess "Error: Unable to open output file: $lgtseek->{output_dir}/$output_prefix\-m8.txt\n";

    my $last_id = 'foo';
    my %bho;
    $bho{$last_id} = 'foobar';

    while (<IN>) {
        my @fields = split(/\t/);

        my $evalue = $fields[10];
        next if ( $evalue >= $evalue_cutoff );

        my $id = $fields[0];
        if ( $fields[0] =~ /^(.*)(\/\d)*/ ) { $id = $1; }

        if ( $bho == 1 ) {
            if ( $bho{$id} ) {
                if ( $evalue eq $bho{$id} ) {
                    $last_id = $id;
                    print OUT;
                }
            }
            elsif ( !$bho{$id} ) {
                delete $bho{$last_id};
                $last_id = $id;
                $bho{$id} = $evalue;
                print OUT;
            }
        }
        else {
            print OUT;
        }
    }
    close IN;
    close OUT;
}
else {
    my $in = Bio::SearchIO->new(
        -file   => "blastall -p blastn -a $lgtseek->{threads} -e $evalue_cutoff -T F -d $options{db} -i $fasta |",
        -format => 'blast',
        -best   => '$bho',
        -signif => $evalue_cutoff,
    );

    my $out = Bio::SearchIO->new(
        -file          => ">$lgtseek->{output_dir}/$output_prefix\_blastn.txt",
        -output_format => 'blast',
    );

    while ( my $result = $in->next_result ) {
        $out->write_result($result);
    }
}
$lgtseek->time_check;
print_complete( \%options );

sub help {
    die "Help: This script will run blast-m8 and filter for best hits &| evalue.
    ----------------------------------------------------------------------------------------------------
        --fasta|f=          <FASTA> File to blastn against the db. 
        --bam|b=             <BAM>  File to blastn against the db. 
    ----------------------------------------------------------------------------------------------------
        --db|d=             <FastaDB> Fasta reference file to blastn against. 
        --formatdb=         <0|1> [0] 1= Build the blast database to blast again. 
    ----------------------------------------------------------------------------------------------------
        --evalue_cutoff =       [1] Max Evalue allowed. Example: 1e-5
        --best_hits_only|bho =      <0|1> [0] 1=Only keep best hits.
        --m8=               <0|1> [1] 1= Output blast -m8 format. 0= \"Normal\" blast output.
    ----------------------------------------------------------------------------------------------------
        --output_dir|o=         Directory for all output. Will be created if it doesn't exist.
          --output_prefix|op=       Name of prefix for the output. 
          --subdirs=            <0|1> [0] 1 = Create a new subdirectory for each input in the output directory.
    ----------------------------------------------------------------------------------------------------
        --Qsub|Q=           <0|1> [0] 1 = Submit the job to the grid.
          --sub_name=           Name of the job to submit.
          --project=            [jdhotopp-lab] Grid project to use.
    ----------------------------------------------------------------------------------------------------    
        --conf_file=            [/home/USER/.lgtseek.conf]
    ----------------------------------------------------------------------------------------------------    
        --verbose|V         <0|1> [1] 1= Verbose reporting of progress. 0 =Turns off reports.
        --print_hostname|ph=        <0|1> [1] Print hostname. Defaults to 1 if verbose=1. 
        --help
    ----------------------------------------------------------------------------------------------------\n";
}

__END__
