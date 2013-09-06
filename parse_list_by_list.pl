#!/usr/bin/perl
use warnings;
use strict;
use lib qw(/local/projects-t3/HLGT/scripts/lgtseek/lib/ /opt/lgtseek/lib/);      ### May need to change this depending on where the script is being run
use LGTSeek;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
		'input_list=s',
		'good_list=s',
		'bad_list=s',
		'output_dir=s',
		'output_prefix=s',
		'output_suffix=s',
		'output_name=s',
		'help',
		);

if($options{help}){die "Help: This script will parse out the good or bad lines(ids etc) from --input_list.
		--input_list=
		--good_list=
		--bad_list=
		--output_dir=
		--output_prefix=
		--output_suffix=
		--output_name=
		--help.\n"
		};

my $lgtseek = LGTSeek->new2({
	options => \%options,
	});

if(!$lgtseek->{input_list}){die "Must give a list to parse for desired reads. use --input_list=<LIST>\n";}
if(!$lgtseek->{good_list} && !$lgtseek->{bad_list}){die "Must pass a --good_list or --bad_list. Try again\n";}
if($lgtseek->empty_chk({input=>$lgtseek->{input_list}})==1){die "Error: --input_list is empty\n";}
my $good_ids = {};
my $bad_ids = {};
if($lgtseek->{good_list}){$good_ids = $lgtseek->_read_ids({ list => $lgtseek->{good_list} });} 
if($lgtseek->{bad_list}){$bad_ids = $lgtseek->_read_ids({ list => $lgtseek->{bad_list} });} 
## Setup input and output bams
my $input = $lgtseek->{input_list};
my ($fn,$path,$suf)=fileparse($input,('/\.\w+$/',".list"));
my $out_dir = $lgtseek->{output_dir} ? $lgtseek->{output_dir} : $path;
my $prefix = $lgtseek->{output_prefix} ? $lgtseek->{output_prefix} : $fn;
my $suffix = $lgtseek->{output_suffix} ? $lgtseek->{output_prefix} : "_filtered.list";
my $out_fn = $lgtseek->{output_name} ? $lgtseek->{output_name} : "$prefix$suffix";
my $out = "$out_dir$out_fn";
open(my $ifh, "<","$lgtseek->{input_list}") or die "Can't open input_list: $lgtseek->{input_list} because: $!\n";
open(my $ofh, ">", "$out") or die "Can't open output: $out because: $out\n";
while(<$ifh>){
	chomp;
	if( $lgtseek->{good_list} && $good_ids->{$_}){print $ofh "$_\n";}
	if( $lgtseek->{bad_list} && !$bad_ids->{$_}){print $ofh "$_\n";}
}
close $ifh;
close $ofh;
print STDERR "Completed parsing: $lgtseek->{input_list} for lines from:$lgtseek->{good_list} $lgtseek->{bad_list}.\n";
