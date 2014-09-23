#!/usr/local/bin/perl

use strict;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Data::Dumper;
use Bio::Index::Fastq;
use IO::String;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options = ();
my $results = GetOptions( \%options, 'id_file=s', 'input_fasta=s', 'format=s', 'index_file=s', 'hash', 'help|h' ) || pod2usage();

# display documentation
if ( $options{'help'} ) {
    die "Help. This script will put out sequences from multiFasta ID.
      --id_file=
      --input_fasta=
      --format=
      --index_file=
      --hash           Use this to pull the id's into a hash and print as reading through the fasta. Skip index fasta file. Also prints both reads if ID overlaps. w/o --hash only prints 2nd read if non-unqiue id. 
      --help\n";
}
########################  SLOW    ############################
## This should only be used with a smaller list/fasta file. ##
##############################################################
if ( $options{hash} ) {
    open( ID, "<", "$options{id_file}" ) or die "Could not open id_file: $options{id_file} because: $!\n";
    my %id;
    while (<ID>) {
        chomp;
        $_ =~ s/^\>//;
        $id{$_}++;
    }
    close ID;
    my $reader = new Bio::SeqIO( -format => 'fasta', -file => $options{input_fasta} );
    while ( my $seqRec = $reader->next_seq ) {
        my $header = $seqRec->display_id();
        $header =~ s/^\>//;
        if ( $id{$header} ) {
            print "\>" . $seqRec->display_id() . "\n" . $seqRec->seq() . "\n";
        }
    }
    exit;
}

######################  FAST   ###############################
########## INDEX Fasta file and pull seq's by id #############
##############################################################

my $index_file = "$options{input_fasta}.idx";
if ( $options{index_file} ) {
    $index_file = $options{index_file};
}
my $make_index = !-e "$index_file";
my $format = $options{format} ? $options{format} : 'fasta';

# create database from directory of fasta files
my $db;
if ( $format eq 'fasta' ) {
    $db = Bio::DB::Fasta->new( $options{input_fasta} );
}

elsif ( $format eq 'fastq' ) {
    my $write = 1;
    if ( -e $index_file ) {
        $write = 0;
    }
    $db = Bio::Index::Fastq2->new(
        '-filename'   => "$index_file",
        '-write_flag' => $write
    );
    if ($make_index) {
        print STDERR "Didn't find the index so we're making it\n";
        $db->make_index( $options{input_fasta} );
    }
}

open IN, "<$options{id_file}" or die;
while (<IN>) {
    chomp;
    next if ( !$_ );

    my $string;
    my $stringio = IO::String->new($string);

    my $out;
    if ( $format eq 'fasta' ) {
        $out = Bio::SeqIO->new(
            '-format' => 'Fasta',
            '-fh'     => $stringio
        );
    }
    elsif ( $format eq 'fastq' ) {
        $out = Bio::SeqIO->new(
            '-format' => 'Fastq',
            '-fh'     => $stringio
        );
    }

    my $seq = $db->get_Seq_by_id($_);

    if ( !$seq ) {
        print STDERR "Couldn't find sequence for $_\n";
        next;
    }
    my $header = $seq->display_id();

    #    my $header = $db->header($_);
    $out->write_seq($seq);

    #    $string =~ s/\>.*/>$header/;
    print $string;

    #print "\n";
}

