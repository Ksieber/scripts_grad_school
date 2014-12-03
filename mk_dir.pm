package mk_dir;
use warnings;
use strict;
use Carp;
$Carp::MaxArgLen = 0;    ## Report full length error
use run_cmd;
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw( mk_dir );

=head2 mk_dir
    Title       : mk_dir
    Usage       : mk_dir($path);
    Function    : Creates the directory with proper permissions
=cut 

sub mk_dir {
    my $dir = shift;
    chomp($dir);
    $dir =~ s/\/{2,}/\//g;
    run_cmd("mkdir -p -m u=rwx,g=rwx,o= $dir");
}
1;
