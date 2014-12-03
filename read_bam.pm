package read_bam;
use warnings;
use strict;
use Carp;
use Exporter;
use parse_flag;
our @ISA    = qw(Exporter);
our @EXPORT = qw( open_bam read_bam );

sub open_bam {
    my $bam = shift;
    if ( !$bam )    { confess "Error: read_bam.pm &open didn't recience a bam properly.\n"; }
    if ( !-e $bam ) { confess "Error: read_bam.pm &open bam doesn't exist.\n"; }
    open( my $FH, "-|", "samtools view $bam" ) or confess "Error: read_bam.pm &open couln't open input: $bam\n";
    return $FH;
}

sub read_bam {
    my $FH   = shift;
    my $line = <$FH>;
    if ( !$line ) { close $FH; return undef; }
    chomp($line);
    my ( $id, $flag, $chr, $position, $mapq, $cigar, $mate_ref, $mate_position, $insert, $sequence, $qual ) = ( split /\t/, $line )[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 ];
    $id =~ s/(.+)\/\d+/$1/;
    my $converted_flag = parse_flag($flag);
    return {
        id            => $id,
        flag          => $converted_flag,
        chr           => $chr,
        position      => $position,
        mapq          => $mapq,
        cigar         => $cigar,
        mate_ref      => $mate_ref,
        mate_position => $mate_position,
        insert        => $insert,
        sequence      => $sequence,
        quality       => $qual,
        line          => $line,
    };
}

1;
