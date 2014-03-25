package parse_flag;
use warnings;
use strict;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( parse_flag );



sub parse_flag {
	my $int = shift;
	my $rawbin = &_dec2bin($int);
	my $rev = scalar $rawbin;
	my $bin = sprintf("%011d", $rev);
	my $final_bin = reverse $bin;
	return {
		'paired' => substr($final_bin, 0, 1),
		'proper' => substr($final_bin, 1, 1),
		'qunmapped' => substr($final_bin, 2, 1),
		'munmapped' => substr($final_bin, 3, 1),
		'qrev' => substr($final_bin, 4, 1),
		'mrev' => substr($final_bin, 5, 1),
		'first' => substr($final_bin, 6, 1),
		'last' => substr($final_bin, 7, 1),
		'secondary' => substr($final_bin, 8, 1),
		'failqual' => substr($final_bin, 9, 1),
		'pcrdup' => substr($final_bin, 10, 1)
	};
}

## Private subroutine
sub _dec2bin {
	my $str = unpack("B32", pack("N", shift));
	$str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
	return $str;
}

1;