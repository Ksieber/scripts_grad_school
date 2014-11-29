#!/usr/bin/perl
use lib ( '/home/ksieber/perl5/lib/perl5/', '/home/ksieber/scripts/' );
use strict;
use warnings;
use lca2krona;
use run_cmd;
use mk_dir;
use setup_input;
use print_call;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'input_list=s', 'output_dir|o=s', 'output_prefix|p=s', 'col=i', 'Qsub|q=i', 'sub_name=s', 'project=s', 'help|?' )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} )   { &help; }
if ( !$options{input} ) { die "Error: You must pass an input file with --input=<LCA.txt>\n"; }

if ( defined $options{Qsub} and $options{Qsub} == 1 ) {
    if ( !$options{sub_name} ) { $options{sub_name} = "lca2krona" }
    if ( !$options{project} )  { $options{project}  = $lgtseek->{project}; }
    Qsub_script( \%options );
}

my $inputs = setup_input( \%options );
foreach my $input (@$inputs) {
    my $krona = lca2krona(
        {   input         => $input,
            output_dir    => $options{output_dir},
            output_prefix => $options{output_prefix},
            col           => $options{col},
        }
    );
}

print_complete( \%options );

sub help {
    die "This script will create an krona.html file from a LCA.txt.
    --input=        </full/path/to/LCA.txt>
    --input_list=       </full/path/to/list_of_LCA.list> 1 file / line.      
    --col=      <#> [1] 1= Use the 2nd tab spaced column to parse the LCAs (zero based).
    --output_dir=       </path/for/output/>     [input file's dir] 
    --output_prefix=    <\$prefix_krona.html>   [input file's name]\n";
}
