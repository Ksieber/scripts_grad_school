package read_in_list;
use strict;
use warnings;
use Carp;
$Carp::MaxArgLen = 0;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( read_in_list );

=head2 &read_in_list

 Title   : read_in_list
 Usage   : 
 Function: 
 Returns : 
 Args    : 
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