package linecount;
use run_cmd;
use strict;
use warnings;
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw( wc );

sub wc {
    my $file = shift;
    my $count;
    if ( $file =~ /.*\.bam$/ ) {
        $count = `samtools view $file | wc -l`;
    }
    elsif ( $file =~ /.*\.sam$/ ) {
        $count = `samtools view -S $file | wc -l`;
    }
    elsif ( $file =~ /.*\.gz$/ ) {
        $count = `zcat $file | wc -l`;
    }
    else {
        $count = `wc -l $file`;
    }
    $count =~ /^(\d+)\s+/;
    my $num = $1;
    if ( $num == 0 ) {
        print STDERR "* Warning * | This file is empty: $file\n";
    }
    return "$num";
}
1;

__END__