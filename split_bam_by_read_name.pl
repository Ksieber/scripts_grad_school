#!/usr/local/bin/perl -w
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options;
my $results = GetOptions( \%options, 'input_bam=s', 'split_name=s', 'output_dir=s', 'dedup', 'split_size=s', 'help|h', );

if ( $options{help} ) {
    die "\nHELP: This script will take large .bam(s) file
and split it into smaller bam file(s).  
    *--input_bam=*        Single .bam file to split. REQUIRED.
    *--split_name=*       Split a single read group out, or a list of COMMA seperated read groups.  REQUIRED
    --output_dir=         [Current Working Directory] unless specified otherwise
    --dedup               Remove flags marked as duplicates (PCR/optical)
    --split_size=         [Infinite] Specify the number of lines per smaller .bam files
    --help\n";
}

if ( !$options{input_bam} && !$options{list_bams} ) { die "ERROR:  Must give an input bam(s) to  apart with --input_bam or --list_bams\n"; }
if ( !$options{split_name} ) { die "ERROR:  Must give an read group name to split on with --split_names\n"; }

my @read_groups = split( /,/, $options{split_name} );
my %output_hash;
my %fh_hash;
my $n = 0;
$options{input_bam} =~ /^(.*)\.bam/;

foreach my $read_grp (@read_groups) {
    $output_hash{$read_grp} = $options{output_dir} ? "$options{output_dir}$1_$read_grp.sam" : "$1_$read_grp.sam";
    open( my $fh, ">", "$output_hash{$read_grp}" ) || die "ERROR: Couldn't open $output_hash{$read_grp} because : $!\n";
    $fh_hash{$read_grp} = $fh;
    my $get_header_command = "samtools view -H $options{input_bam}";
    open( my $grab_header, "$get_header_command |" );
    while (<$grab_header>) {
        print { $fh_hash{$read_grp} } $_;
    }
}

my $samtools_command = $options{dedup} ? "samtools view -F 0x400 $options{input_bam}" : "samtools view $options{input_bam}";
open my $samtools_view, "$samtools_command |";
while (<$samtools_view>) {
    chomp;
    my @fields = split;
    foreach my $read_grp (@read_groups) {
        if ( $fields[0] =~ /$read_grp/ ) {
            my $output_fields = join( "\t", @fields[ 0 .. 10 ] );
            print { $fh_hash{$read_grp} } "$output_fields\n";
        }
    }
}

foreach my $read_grp (@read_groups) {
    close $read_grp;
    $output_hash{$read_grp} =~ /^(.*)\.sam/;
    my $output_bam = "$1.bam";
    `samtools view -bS $output_hash{$read_grp} > $output_bam`;
    `rm $output_hash{$read_grp}`;
}

