package print_call;
use warnings;
use strict;
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw( print_call print_complete print_hostname print_notebook print_progress);

=head2 print_call
    Title       : print_call
    Usage       : print_call(\%options,"other comments", "to also print", "like version number");
    Function    : Print STDERR original perl script & --options with it
    Args        : Hash ref
=cut 

sub print_call {
    my $options = shift;
    print STDERR "================================================================================================================================\n";
    print STDERR "$0";
    foreach my $key ( sort keys %$options ) {
        if ( $options->{$key} ) { print STDERR " --$key=$options->{$key}" }
    }
    if ( @_ and scalar(@_) >= 1 ) {
        foreach my $comment (@_) {
            print STDERR "\n$comment";
        }
    }
    print STDERR "\n================================================================================================================================\n";
}

sub print_complete {
    my $options = shift;
    print STDERR "================================================================================================================================\n";
    print STDERR "$0";
    foreach my $key ( sort keys %$options ) {
        if ( $options->{$key} ) { print STDERR " --$key=$options->{$key}" }
    }
    if ( @_ and scalar(@_) >= 1 ) {
        foreach my $comment (@_) {
            print STDERR "\n$comment";
        }
    }
    print STDERR "\n--------------------------------------------------------------------------------------------------------------------------------\n";
    print STDERR "++++++++++ COMPLETED:\t$0\t:++++++++++\n";
    print STDERR "================================================================================================================================\n";
}

sub print_hostname {
    my $options = shift;
    if ( $options->{print_hostname} ) {
        print STDERR "================================================================================================================================\n";
        system("echo =================       HOSTNAME: `hostname` 1>&2");
    }
}

sub print_notebook {
    my $options = shift;
    if ( $options->{output_dir} ) {
        open( OUT, ">>", "$options->{output_dir}/\_notebook.txt" )
            || die "Error: Unable to open output: $options->{output_dir}\_notebook.txt because: $!\n";
    }
    else {
        open( OUT, ">", "_notebook.txt" ) || die "Error: Unable to open output: _notebook.txt because: $!\n";
    }
    print OUT "$0";
    foreach my $key ( sort keys %$options ) {
        if ( $options->{$key} ) { print OUT " --$key=$options->{$key}" }
    }
    if ( @_ and scalar(@_) >= 1 ) {
        foreach my $comment (@_) {
            print OUT "\n$comment";
        }
    }
    print OUT "\n";
    close OUT;
}

=head2 print_progress
    Title       : print_progress
    Usage       : while(<>){print_progress; do important stuff ... }
    Function    : Print "." to STDERR every $x_#_of_iterations through a loop to show progress is being made.
    Args        : print_progress( $number )
                    default = 100000. ie every 100000th iteration through a loop, "." is printed.
=cut

my $counter;

sub print_progress {
    my $cut_off_opt = shift;
    my $cut_off = defined $cut_off_opt ? $cut_off_opt : 100000;
    $counter++;
    if ( $counter >= ( $cut_off * 10 ) ) { $counter = 0; print STDERR "\r                                                                                        "; }
    elsif ( $counter % $cut_off == 0 ) {
        my $string;
        for ( my $n = 1; $n <= ( $counter / $cut_off ); $n++ ) { $string = $string . ".\t"; }
        print STDERR "\r$string";
    }
}

1;
