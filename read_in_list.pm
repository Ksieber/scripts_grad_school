package read_in_list;
use strict;
use warnings;
use Carp;
$Carp::MaxArgLen = 0;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( read_in_list );

=head2 &empty_chk

 Title   : empty_chk
 Usage   : if(empty_chk($file)== 1){die "File is empty\n";}
 Function: Check to make sure a file is not empty
 Returns : 1 = Empty; 0 = input is not empty
 Args    : File to check;
=cut
sub read_in_list {
    my $file = shift;
    if(!$file){confess "Must pass &read_in_list a file."}
    my @ret;
    open(IN,"<","$file") or confess "Error: Not able to open the list file: $file\n";
    chomp(@ret = <IN>);
    close IN;
    return \@ret;   
}