package setup_input;
use strict;
use warnings;
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw( setup_input );
### Returns array(ref) of input files
## Suggested use: my $input=setup_input(\%options);
## Relies on main script using $options{input} &| $options{input_list}

sub setup_input {
    my @inputList;
    my $options = shift;
    if ( !$options->{input} && !$options->{input_list} ) {
        die "Error: setup_input.pm requires either \$options{input} or \$options{input_list}.\n";
    }
    if ( $options->{input} ) {
        push( @inputList, $options->{input} );
    }
    if ( $options->{input_list} ) {
        open( LIST, "<", "$options->{input_list}" )
            or die "Can't open --input_list: $options->{input_list} because: $!\n";
        while (<LIST>) {
            chomp;
            push( @inputList, $_ );
        }
        close LIST;
    }
    return \@inputList;
}
1;
