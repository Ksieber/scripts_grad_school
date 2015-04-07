#!/usr/local/bin/perl
use warnings;
use strict;
use File::Basename;
use Bio::SeqIO;

if ( -t STDIN && !@ARGV ) { &help; }

## Global Variables
my %counts;
&clear_counts;
my $previous_id;

## Determine how to process the input
if ( defined $ARGV[0] ) {
    if ( -e $ARGV[0] ) {
        my ( $fn, $dir, $suf ) = fileparse( $ARGV[0], ( '.fa', '.fasta', '.fq', '.fastq', '.bam' ) );
        if ( $suf eq ".fa" || $suf eq ".fasta" ) { &process_fasta( { file => $ARGV[0] } ); }
        elsif ( $suf eq ".fq" || $suf eq ".fastq" ) { &process_fastq( { file => $ARGV[0] } ); }
        elsif ( $suf eq ".bam" ) { &process_bam( { file => $ARGV[0] } ); }
    }
    elsif ( $ARGV[0] =~ /[AaTtGgCc]+/ ) { &count( $ARGV[0] ); &print_counts; }
}
else {
    chomp( my $line = <STDIN> );
    my @split_line = split( /\t/, $line );
    if ( $line =~ /^\>/ ) { &process_fasta( { stdin => 1 } ); }
    elsif ( $line =~ /^\@/ ) { &process_fastq( { stdin => 1 } ); }
    elsif ( @split_line >= 10 ) { &process_bam( { line => $line } ); }
    elsif ( $line =~ /[AaTtGgCc]+/ ) { &count($line); &print_counts; }
    else                             { die "Could not determine what type of file this is. Please try again.\n"; }
}

## Subroutines
sub process_bam {
    my $opts = shift;
    if ( defined $opts->{file} ) {
        open( IN, "-|", "samtools view $opts->{file}" ) or die "Unable to open the input file: $opts->{file}\n";
        while (<IN>) {
            chomp( my $bam_line = $_ );
            &_process_bam_line($bam_line);
        }
        close IN;
        &print_grand_total;
    }
    elsif ( defined $opts->{line} ) {
        my $bam_line = $opts->{line};
        &_process_bam_line($bam_line);
    }
    else {
        die "Unable to process_bam properly.\n";
    }
}

sub _process_bam_line {
    my $bam_line = shift;
    my ( $id, $sequence ) = ( split /\t/, $bam_line )[ 0, 9 ];
    &count($sequence);
    print "$id: ";
    &print_counts;
    &clear_counts;
}

sub count {
    my $sequence = shift;
    $counts{A}           += $sequence =~ tr/Aa//;
    $counts{T}           += $sequence =~ tr/Tt//;
    $counts{G}           += $sequence =~ tr/Gg//;
    $counts{C}           += $sequence =~ tr/Cc//;
    $counts{N}           += $sequence =~ tr/Nn//;
    $counts{total}       += $sequence =~ tr/AaTtGgCcNn//;
    $counts{tot_A}       += $counts{A};
    $counts{tot_T}       += $counts{T};
    $counts{tot_G}       += $counts{G};
    $counts{tot_C}       += $counts{C};
    $counts{tot_N}       += $sequence =~ tr/Nn//;
    $counts{grand_total} += $counts{total};
}

sub clear_counts {
    $counts{A}     = 0;
    $counts{T}     = 0;
    $counts{G}     = 0;
    $counts{C}     = 0;
    $counts{N}     = 0;
    $counts{total} = 0;
}

sub print_counts {
    print "\%A: ";
    printf "%5.1f\t", ( ( $counts{A} / $counts{total} ) * 100 );
    print "\%T: ";
    printf "%5.1f\t", ( ( $counts{T} / $counts{total} ) * 100 );
    print "\%G: ";
    printf "%5.1f\t", ( ( $counts{G} / $counts{total} ) * 100 );
    print "\%C: ";
    printf "%5.1f\t", ( ( $counts{C} / $counts{total} ) * 100 );
    if ( $counts{N} >= 1 ) { print "\%N: "; printf "%2.1f\t", ( ( $counts{N} / $counts{total} ) * 100 ); }
    print "\%GC: ";
    printf "%5.1f\t", ( ( ( $counts{G} + $counts{C} ) / $counts{total} ) * 100 );
    print "%AG: ";
    printf "%5.1f\t", ( ( ( $counts{A} + $counts{G} ) / $counts{total} ) * 100 );
    print "Total_bp: $counts{total}\n";
}

sub print_grand_total {
    print "_____________________________________________________________________________________________________________\n";
    print "Grand_Total: ";
    print "\%A: ";
    printf "%2.1f\t", ( $counts{tot_A} / $counts{grand_total} ) * 100;
    print "\%T: ";
    printf "%2.1f\t", ( $counts{tot_T} / $counts{grand_total} ) * 100;
    print "\%G ";
    printf "%2.1f\t", ( $counts{tot_G} / $counts{grand_total} ) * 100;
    print "\%C ";
    printf "%2.1f\t", ( $counts{tot_C} / $counts{grand_total} ) * 100;
    if ( $counts{tot_N} >= 1 ) { print "\%N:"; printf "%2.1f", ( $counts{tot_N} / $counts{grand_total} ) * 100; }
    print "\%GC: ";
    printf "%2.1f\t", ( ( $counts{tot_G} + $counts{tot_C} ) / $counts{grand_total} ) * 100;
    print "%AG: ";
    printf "%2.1f\t", ( ( $counts{tot_A} + $counts{tot_G} ) / $counts{grand_total} ) * 100;
    print "Grand_Total_bp: $counts{grand_total}\n";
    print "_____________________________________________________________________________________________________________\n";
}

sub process_fasta {
    my $opts = shift;
    my $seq_fh;
    if ( defined $opts->{file} ) {
        $seq_fh = Bio::SeqIO->new( -format => 'Fasta', -file => $opts->{file} );
    }
    elsif ( $opts->{stdin} == 1 ) {
        $seq_fh = Bio::SeqIO->new( -format => 'Fasta', -fh => \*STDIN );
    }
    else {
        die "Unable to read the fasta properly.\n";
    }
    &process_seq_fh($seq_fh);
}

sub process_fastq {
    my $opts = shift;
    my $seq_fh;
    if ( defined $opts->{file} ) {
        $seq_fh = Bio::SeqIO->new( -format => 'Fastq', -file => $opts->{file} );
    }
    elsif ( $opts->{stdin} == 1 ) {
        $seq_fh = Bio::SeqIO->new( -format => 'Fastq', -fh => \*STDIN );
    }
    else {
        die "Unable to read the fasta properly.\n";
    }
    &process_seq_fh($seq_fh);
}

sub process_seq_fh {
    my $seq_fh    = shift;
    my $seqs_seen = 0;
    while ( my $seq = $seq_fh->next_seq() ) {
        print $seq->display_id() . ": ";
        count( $seq->seq() );
        &print_counts;
        &clear_counts;
        $seqs_seen++;
    }
    &print_grand_total if ( $seqs_seen > 1 );
}

sub help {
    die "
    _____________________________________________________________________________________________
    This script will calculate the percentage of each residue, GC, and AG of a sequence.
    It can read sequence, fasta files, and bam files from STDIN or \$ARGV[0].
    _____________________________________________________________________________________________\n";
}
