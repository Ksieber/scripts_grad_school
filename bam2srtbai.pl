#!/usr/local/bin/perl
use warnings;
use strict;
use File::Basename;
use errorcheck;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
      'input=s',
      'output_prefix=s',
      'output_dir',
      'help',
);

if($options{help}){
   die "Help: This script will take a bam to sort and index it.
      --input=           <BAM> to sort & index
      --output_prefix=  prefix.bam
      --output_dir=     directory to put output
      --help\n";
}

my $input = $options{input} ? $options{input} : $ARGV[0];
my ($bam,$path,$suf)=fileparse($input,".bam");
my $prefix = $options{output_prefix} ? $options{output_prefix} : $bam;
my $dir = $options{output_dir} ? $options{output_dir} : $path;
$dir =~ s/\/$//;
my $out = "$dir/$prefix";

`samtools sort $input $out\.srt`;
&errchk($?);
`samtools index $out\.srt.bam $out\.srt.bai`;
&errchk($?);


