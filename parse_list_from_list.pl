#!/usr/bin/perl
use warnings;
no warnings 'uninitialized';
use strict;
use lib qw(/local/projects-t3/HLGT/scripts/lgtseek/lib/ /opt/lgtseek/lib/);    ### May need to change this depending on where the script is being run
use LGTSeek;
use File::Basename;
use Cwd;

if ( !@ARGV ) { &help; }

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input_list|I=s', 'input|i=s', 'good_list=s', 'bad_list=s', 'output_dir|o=s', 'output_prefix|p=s', 'output_suffix|s=s', 'output_name=s', 'help|?', );

if ( $options{help} ) { &help; }

my $lgtseek = LGTSeek->new2( \%options );

if ( !$lgtseek->{input_list} and !$lgtseek->{input} ) { die "Must give a list to parse for desired reads. use --input_list|I=<LIST>\n"; }
if ( !$lgtseek->{good_list} && !$lgtseek->{bad_list} ) { die "Must pass a --good_list or --bad_list. Try again\n"; }
my $good_ids = {};
my $bad_ids  = {};
if ( $lgtseek->{good_list} ) { $good_ids = $lgtseek->_read_ids( { list => $lgtseek->{good_list} } ); }
if ( $lgtseek->{bad_list} )  { $bad_ids  = $lgtseek->_read_ids( { list => $lgtseek->{bad_list} } ); }
## Setup input and output bams
my $input = defined $lgtseek->{input_list} ? $lgtseek->{input_list} : $lgtseek->{input};
if ( $lgtseek->empty_chk( { input => $input } ) == 1 ) { die "Error: --input_list is empty\n"; }
my ( $fn, $path, $suf ) = fileparse( $input, ( '.list', '.txt', '\.\w+$' ) );
my $out_dir = $lgtseek->{output_dir}    ? $lgtseek->{output_dir}    : getcwd;
my $prefix  = $lgtseek->{output_prefix} ? $lgtseek->{output_prefix} : $fn;
my $suffix  = $lgtseek->{output_suffix} ? $lgtseek->{output_suffix} : "_filtered.list";
my $out_fn  = $lgtseek->{output_name}   ? $lgtseek->{output_name}   : "$prefix$suffix";
my $out     = "$out_dir$out_fn";

open( my $ifh, "<", "$input" ) or die "Can't open input_list: $input because: $!\n";
my $ofh;
if ( defined $options{output_dir} or defined $options{output_prefix} or defined $options{output_suffix} or defined $options{output_name} ) {
    open( $ofh, ">", "$out" ) or die "Can't open output: $out because: $out\n";
}
else {
    $ofh = *STDOUT;
}

while (<$ifh>) {
    chomp;
    if ( $lgtseek->{good_list} && $good_ids->{$_} ) { print $ofh "$_\n"; }
    if ( $lgtseek->{bad_list}  && !$bad_ids->{$_} ) { print $ofh "$_\n"; }
}
close $ifh;
close $ofh;
print STDERR "\n\tCompleted parsing: $input for lines from:" . "$lgtseek->{good_list}" . "$lgtseek->{bad_list}" . ".\n\n";

sub help {
    die "\nHelp: This script will parse the input for the good or bad lines(ids etc) from --input_list.
		--input_list|I=
		--good_list=
		--bad_list=
		--output_dir|o=
		--output_prefix|p=
		--output_suffix|s=
		--output_name=
		# If no --output given, output goes to STDOUT\n\n"
}
