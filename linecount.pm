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
    elsif ( $file =~ /.*\.gz$/ ) {
        $count = `zcat $file | wc -l`;
    }
    else {
        $count = `wc -l $file`;
    }
    $count =~ /^(\d+)\s+/;
    my $num = $1;
    if ( $num == 0 ) {
        print STDERR "Warning: The file is empty. File:$file\tLines:$num\n";
    }
    return "$num";
}
1;

__END__

# sub wc {
#    my $file = shift;
#    my $count=`wc -l $file`;
#    $count=~/^(\d+)\s+/;
#    my $num=$1;
#    if($num==0){
#       print STDERR "Warning: The file is empty. File:$file\tLines:$num\n";
#    }
#    return "$num";
# }
