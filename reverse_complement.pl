#!/usr/local/bin/perl
use warnings;
use strict;
use Bio::Util::DNA qw(:all);
use Bio::SeqIO;
use File::Basename;

if ( -t STDIN and not @ARGV ) { die "\nPlease pass a sequence to this script through STDIN or ARGV.\nExample: reverse_complement.pl ATGCATGC\n\n"; }

if ( defined $ARGV[0] ) {
    if ( -e $ARGV[0] ) {
        my ( $fn, $dir, $suf ) = fileparse( $ARGV[0], ( '.fa', '.fasta', '.fq', '.fastq' ) );
        if ( $suf eq ".fa" || $suf eq ".fasta" ) { &process_fasta( { file => $ARGV[0] } ); }
        elsif ( $suf eq ".fq" || $suf eq ".fastq" ) { &process_fastq( { file => $ARGV[0] } ); }
    }
    elsif ( $ARGV[0] =~ /[AaTtGgCc]+/ ) { &process_seq_text( $ARGV[0] ); }
}
else {
    chomp( my $line = <STDIN> );
    if ( $line =~ /^\>/ ) { &process_fasta( { stdin => 1 } ); }
    elsif ( $line =~ /^\@/ ) { &process_fastq( { stdin => 1 } ); }
    elsif ( $line =~ /[AaTtGgCc]+/ ) { &process_seq_text($line); }
    else                             { die "Could not determine what type of file this is. Please try again.\n"; }
}

## Subroutines
sub process_seq_text {
    my $sequence = shift;
    my $rev_comp = reverse_complement( \$sequence );
    print "$$rev_comp\n";
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
    my $seq_fh = shift;
    while ( my $seq = $seq_fh->next_seq() ) {
        printf "%-40s : %-100s\n", ( $seq->display_id(), $seq->revcom()->seq() );
    }
}

