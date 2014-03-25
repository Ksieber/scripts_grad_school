=head1 SYNOPSIS

	This module allows printing the same message to multiple file handles at once.

	my $mfh = MultipleFileHandles->new([$fh1,$fh2]); ## Can also use ->new(\@list_of_fhs);
	$mfh->add($fh3);
	$mfh->print("Something really important\n");
	$mfh->close;

=cut
package MultipleFileHandles;
use warnings;
use strict;
use Carp;
$Carp::MaxArgLen = 0; ## Report full length error 

sub new {
	my($class,$fh_array_ref) = @_;
	my $fhs = $fh_array_ref;
	bless $fhs;
	return $fhs;
}

sub add {
	my ($fhs,$fh_to_add)=@_;
	push(@{$fhs},$fh_to_add);
	return $fhs;
}

sub print {
	my ($fhs,$msg)=@_;
	foreach my $fh (@{$fhs}){
		print $fh "$msg";
	}
}

sub close {
	my $fhs = shift;
	foreach my $fh (@{$fhs}){
		close $fh or confess "Error: Can't close filehandle because: $!\n";
	}
}

1;













__END__

package Tie::MultiHandle; # do not use directly.
sub TIEHANDLE
{
	my $pkg = shift;
	bless [ @_ ], $pkg
}

sub AUTOLOAD
{
	my( $self, @args ) = @_;
	our( $AUTOLOAD );
	$AUTOLOAD =~ s/.*::(.*)/\L$1/;
	for my $fh ( @$self ) { $fh->$AUTOLOAD( @args ); }
	# $_->$AUTOLOAD(@args) for @$self;
}

sub DESTROY
{

}


package IO::MultiHandle;
use Symbol;

 # takes a list of handles
sub new {
	my $pkg = shift;
	my $self = gensym;
	tie *$self, Tie::MultiHandle::, @_;
	bless $self, $pkg
}