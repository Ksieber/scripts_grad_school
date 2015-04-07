package read_bam;
use warnings;
use strict;
use Carp;
use run_cmd;
use Exporter;
use parse_flag;
our @ISA    = qw(Exporter);
our @EXPORT = qw( open_bam read_bam write_bam );

=head2 open_bam

    Title   :   open_bam
    Usage1   :   my ($header,$fh) = open_bam('/path/to/some.bam','chr:1-100');
    Usage2   :   my $fh = open_bam('/path/to/some.bam');
    Function:   Open filehandle to a bam
    Returns :   (bam-header, bam-filehandle)
    Args    :   First arg MUST be a bam
                Second arg may be a region, must be chr:1-100 format.
=cut

sub open_bam {
    chomp( my $bam = shift );
    my $region = shift;
    if ( !$bam )    { confess "Error: read_bam.pm &open didn't recience a bam properly.\n"; }
    if ( !-e $bam ) { confess "Error: read_bam.pm &open_bam's input: $bam doesn't exist.\n"; }
    my $header = run_cmd("samtools view -H $bam");
    if ( defined $region and $region =~ /[A-Za-z0-9\_\-]+/ ) { $bam = $bam . " \'$region\'"; }
    open( my $FH, "-|", "samtools view $bam" ) or confess "Error: read_bam.pm could not open input: $bam\n";
    return $header, $FH;
}

=head2 read_bam

    Title   :   read_bam
    Usage   :   while(my $read=read_bam($fh) ) {do stuff ...}
    Usage   :   while(my $read=read_bam($fh,0)){do stuff WITHOUT flags parsing ...}
    Function:   Read a bam filehandle, returning an obj with easy access to data and a "parsed flag"
    Function:   Parsing the flags takes a significant amount of time and can be turned off by passing 0 after the filehandle.
    Args    :   Filehandle, (<0|1>[1])
    Returns :   An object with keys (below) to access the bam data: $object->{$key}=bam_data_value
                    id
                    flag->{ paired proper qunmapped munmapped qrev mrev first last secondary failqual pcrdup supplementary }
                    chr
                    position
                    mapq
                    cigar
                    mate_ref
                    mate_position
                    insert
                    sequence
                    quality
                    line
 
 TODO  :   Add "XA" tag access

=cut

sub read_bam {
    my $FH         = shift;
    my $parse_flag = shift;
    my $line       = <$FH>;
    if ( !$line ) { close $FH; return undef; }
    chomp($line);
    my ( $id, $flag, $chr, $position, $mapq, $cigar, $mate_ref, $mate_position, $insert, $sequence, $qual ) = ( split /\t/, $line )[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 ];
    $id =~ s/(.+)\/\d+/$1/;
    my $converted_flag = ( defined $parse_flag and $parse_flag == 0 ) ? $flag : parse_flag($flag);
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

=head2 write_bam

    Title   :   write_bam
    Usage   :   my $out_fh = write_bam('/path/to/some.bam',$header);
    Function:   Open filehandle to a bam
    Returns :   Filehandle
    Args    :   First arg MUST be the full output name.
                Second arg must be the header. This can be the $header returned from open_bam.
=cut

sub write_bam {
    chomp( my $output = shift );
    chomp( my $header = shift );
    open( my $FH, "| samtools view - -bo $output" ) or confess "Error: unable to open output: $output\n";
    print $FH "$header\n";
    return $FH;
}
1;
