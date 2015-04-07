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
 Usage   : my @list = read_in_list($file);
 Function: opens a file, reads in line, and returns an array of the lines
 Returns : ARRAY
 Args    : n/a
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
