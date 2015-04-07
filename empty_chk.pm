package empty_chk;
use strict;
use warnings;
use Carp;
$Carp::MaxArgLen = 0;
use run_cmd;
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw( empty_chk );

=head2 &empty_chk

 Title   : empty_chk
 Usage   : if(empty_chk($file)== 1){die "File is empty\n";}
 Function: Check to make sure a file is not empty
 Returns : 1 = Empty; 0 = input is not empty
 Args    : File to check;
=cut

sub empty_chk {
    my $file = shift;
    if ( !$file ) { confess "Must pass &empty_chk a file." }
    my $empty = 0;    ## 0 = False, 1=True, file is empty.
    my $count;
    if ( $file =~ /\.bam$/ ) {
        $count = run_cmd( "samtools view $file | head | wc -l", 0 );
    }
    elsif ( $file =~ /\.sam$/ ) {
        $count = run_cmd( "samtools view -S $file | head | wc -l", 0 );
    }
    elsif ( $file =~ /\.gz/ ) {
        $count = run_cmd( "zcat $file | head | wc -l", 0 );
    }
    else {
        $count = run_cmd( "head $file | wc -l", 0 );
    }
    if ( $count == 0 ) { $empty = 1; }
    return $empty;
}
