package setup_input;
use strict;
use warnings;
use Exporter;
use read_in_list;
our @ISA    = qw(Exporter);
our @EXPORT = qw( setup_input );
### Returns array(ref) of input files
## Suggested use: my $input=setup_input(\%options);
## Relies on main script using $options{input} &| $options{input_list}

sub setup_input {
    my @inputList;
    my $options = shift;
    if ( !$options->{input} && !$options->{input_list} ) {
        die "Error: &setup_input.pm requires either \$options{input} or \$options{input_list}.\n";
    }
    if ( $options->{input} ) {
        push( @inputList, $options->{input} );
    }
    if ( $options->{input_list} ) {
        if ( !( -e $options->{input_list} ) and $options->{input_list} =~ /\,+/ ) {
            push( @inputList, split( /,/, join( ',', $options->{input_list} ) ) );
        }
        elsif ( -e $options->{input_list} ) {
            push( @inputList, @{ &read_in_list( $options->{input_list} ) } );
        }
    }
    return \@inputList;
}
1;
