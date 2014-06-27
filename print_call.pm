package print_call;
use warnings;
use strict;
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw( print_call print_complete print_hostname print_notebook );

=head2 print_call
	Title		: print_call
	Usage		: print_call(\%options);
	Function	: Print STDERR original perl script & --options with it
	Args		: Hash ref
=cut 

sub print_call {
    my $options = shift;
    print STDERR
        "================================================================================================================================\n";
    print STDERR "$0";
    foreach my $key ( sort keys %$options ) {
        if ( $options->{$key} ) { print STDERR " --$key=$options->{$key}" }
    }
    print STDERR
        "\n================================================================================================================================\n";
}

sub print_complete {
    my $options = shift;
    print STDERR
        "================================================================================================================================\n";
    print STDERR "++++++++++ COMPLETED:\t$0 \t\t++++++++++\n";
    print STDERR
        "--------------------------------------------------------------------------------------------------------------------------------\n";
    print STDERR "$0";
    foreach my $key ( sort keys %$options ) {
        if ( $options->{$key} ) { print STDERR " --$key=$options->{$key}" }
    }
    print STDERR
        "\n================================================================================================================================\n";
}

sub print_hostname {
    my $options = shift;
    if ( $options->{print_hostname} ) {
        print STDERR
            "================================================================================================================================\n";
        system("echo =================       HOSTNAME: `hostname` 1>&2");
    }
}

sub print_notebook {
    my $options = shift;
    if ( $options->{output_dir} ) {
        open( OUT, ">", "$options->{output_dir}/\_notebook.txt" )
            || die "Error: Unable to open output: $options->{output_dir}\_notebook.txt because: $!\n";
    }
    else {
        open( OUT, ">", "_notebook.txt" ) || die "Error: Unable to open output: _notebook.txt because: $!\n";
    }
    print OUT "$0";
    foreach my $key ( sort keys %$options ) {
        if ( $options->{$key} ) { print OUT " --$key=$options->{$key}" }
    }
    print OUT "\n";
    close OUT;
}

1;
