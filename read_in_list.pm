package read_in_list;
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use strict;
use warnings;
use Carp;
$Carp::MaxArgLen = 0;
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw( read_in_list hash_in_list hash_in_data );

=head2 &read_in_list

 Title   : read_in_list
 Usage   : my @list = @{&read_in_list($file)};   ## or my $list = read_in_list($file);
 Function: opens a file, reads in line, and returns an array of the lines
 Returns : ARRAY REF
 Args    : n/a
=cut

sub read_in_list {
    my $file = shift;
    if ( !$file ) { confess "Must pass &read_in_list a file." }
    my @ret;
    open( IN, "<", "$file" ) or confess "Error: Not able to open the list file: $file\n";
    chomp( @ret = <IN> );
    close IN;
    return \@ret;
}

=head2 &hash_in_list

 Title   : hash_in_list
 Usage   : my %hash = %{&read_in_list($file)};   ## or my $hash_ref = read_in_list($file);
 Function: opens a file, reads in line, and returns a hash_ref with keys = line, value = count of exact same lines
 Returns : HASH REF
 Args    : n/a
=cut

sub hash_in_list {
    my $file = shift;
    if ( !$file )    { confess "Must pass &read_in_list a file." }
    if ( !-e $file ) { confess "The input list doesn't exist: $file"; }
    my %ret;
    open( IN, "<", "$file" ) or confess "Error: Not able to open the list file: $file\n";
    while (<IN>) {
        chomp( my $line = $_ );
        $ret{$line}++;
    }
    close IN;
    return \%ret;
}

=head2 &hash_in_data

 Title   : hash_in_data
 Usage   : my %hash = %{&hash_in_data($file)};   ## or my $hash_ref = hash_in_data($file);
 Function: opens a file, reads in line, and returns a hash_ref->{ column0 } = column#.
 Returns : HASH REF
 Args    : Second arguement is optional, specifies the column# to use. Defaults to column1 (0 based counting).  
=cut

sub hash_in_data {
    my $file = shift;
    if ( !$file )    { confess "Must pass &read_in_list a file." }
    if ( !-e $file ) { confess "The input list doesn't exist: $file"; }

    my %ret;
    open( IN, "<", "$file" ) or confess "Error: Not able to open the list file: $file\n";
    while (<IN>) {
        chomp( my $line = $_ );
        my @split_line = split( /\t/, $line );
        if ( defined $split_line[0] ) { $split_line[0] =~ s/(\s+)$//; }
        if ( defined $split_line[1] ) { $split_line[1] =~ s/(\s+)$//; }
        $ret{ $split_line[0] } = $split_line[1];
    }
    close IN;
    return \%ret;
}
