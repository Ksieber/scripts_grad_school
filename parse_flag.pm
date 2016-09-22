package parse_flag;
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use warnings;
use strict;
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw( parse_flag create_flag );

sub parse_flag {
    my $int       = shift;
    my $rawbin    = &_dec2bin($int);
    my $rev       = scalar $rawbin;
    my $bin       = sprintf( "%012d", $rev );
    my $final_bin = reverse $bin;
    return {
        'paired'        => substr( $final_bin, 0,  1 ),
        'proper'        => substr( $final_bin, 1,  1 ),
        'qunmapped'     => substr( $final_bin, 2,  1 ),
        'munmapped'     => substr( $final_bin, 3,  1 ),
        'qrev'          => substr( $final_bin, 4,  1 ),
        'mrev'          => substr( $final_bin, 5,  1 ),
        'first'         => substr( $final_bin, 6,  1 ),
        'last'          => substr( $final_bin, 7,  1 ),
        'secondary'     => substr( $final_bin, 8,  1 ),
        'failqual'      => substr( $final_bin, 9,  1 ),
        'pcrdup'        => substr( $final_bin, 10, 1 ),
        'supplementary' => substr( $final_bin, 11, 1 ),
    };
}

## Private subroutine
sub _dec2bin {
    my $str = unpack( "B32", pack( "N", shift ) );
    $str =~ s/^0+(?=\d)//;    # otherwise you'll get leading zeros
    return $str;
}

sub create_flag {
    my $data = shift;       ## Must be a parsed_flag object (or made to be similar to above data structure lines #15-26)
    my $flag = 0;
    $flag += $data->{paired} *        ( 2**0 );
    $flag += $data->{proper} *        ( 2**1 );
    $flag += $data->{qunmapped} *     ( 2**2 );
    $flag += $data->{munmapped} *     ( 2**3 );
    $flag += $data->{qrev} *          ( 2**4 );
    $flag += $data->{mrev} *          ( 2**5 );
    $flag += $data->{first} *         ( 2**6 );
    $flag += $data->{last} *          ( 2**7 );
    $flag += $data->{secondary} *     ( 2**8 );
    $flag += $data->{failqual} *      ( 2**9 );
    $flag += $data->{pcrdup} *        ( 2**10 );
    $flag += $data->{supplementary} * ( 2**11 );

    return $flag;
}

1;
