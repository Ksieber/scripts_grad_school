#!/usr/local/bin/perl -w
## This script is for spliting up .bam files aligned to human into putative LGT reads and bacterial reads w/ --bac. 

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
			  'input_bam=s',
			  'bac',
			  'output_prefix=s',
			  'output_dir=s',
			  'cleanup',
			  'help|h',
			  );

if ($options{help}){
    die "This script will split up a human aligned .bam file into putative LGT reads.
\t--input_bam\tREQUIRED.  Human aligned.bam to split.
\t--bac\tSplit out UM_UM reads for putative bacterial reads.
\t--output_prefix=\tPrefix for each output.  Ie. SRA_LGT_NC001234
\t--output_dir=\tWill output to current working directory unless another is specified with this. ie. /ksieber_dir/tmp/
\t--cleanup\tWill remove non-essential files after finishing.\n";
}

if (!$options{input_bam}){die "ERROR: No input.  Use --input_bam to specify which human aligned.bam to split.\n";}
if ($options{input_bam}){print STDERR "Starting to process $options{input_bam}\n";}

 my $output_prefix= $options{output_dir} ? "$options{output_dir}/$options{output_prefix}" : "$options{output_prefix}";

## Grab Header for merging files later
print STDERR "Parsing:\n  Header . . .\n";
`samtools view -H $options{input_bam} > $output_prefix.header`;
if ($? != 0 ){die "Fatal Error line 29.\n";}

## Parse out UM_M and M_UM reads; then merge the two files together
print STDERR "  M_UM . . .\n";
`samtools view -b -F0x604 -f0x8 $options{input_bam} > $output_prefix.F4.f8.bam`;
if ($? != 0 ){die "Fatal Error line 31.\n";}
print STDERR "  UM_M . . .\n";
`samtools view -b -F0x608 -f0x4 $options{input_bam} > $output_prefix.F8.f4.bam`;
if ($? != 0 ){die "Fatal Error line 33.\n";}
print STDERR "Merging UM_M and M_UM . . .\n";
`samtools merge -nfh $output_prefix.header $output_prefix.Ff-merged.bam $output_prefix.F4.f8.bam $output_prefix.F8.f4.bam`;
if ($? != 0 ){die "Fatal Error line 36.\n";}

## UM_M and M_UM Are normally not equal causing there to be ~SE reads in Ff-merge.  Remove ~SE reads.
print STDERR "Removing ~SE reads . . .\n";
my $last_read = "dummy";
my $last_line = "big dummy";
open (OUT, ">", "$output_prefix.Ff-merged.PE-filtered.sam") || die "ERROR: Couldn't open $output_prefix.Ff-merged.PE-filtered.sam b/c: $!\n";
my $cmd = "samtools view $output_prefix.Ff-merged.bam | sort -k1";
open (my $PE_filter, "$cmd |");
while (<$PE_filter>){
	chomp;
	if ($_ !~ /^@/){
		my $current_line = $_;
		my $current_read = (split)[0];
		if ($current_read eq $last_read){
			print OUT "$last_line\n$current_line\n";
		}
		if ($current_read ne $last_read){
			$last_read = $current_read;
			$last_line = $current_line;
		}
	}
}
close $PE_filter;
close OUT;

## Convert the PE filtered .sam into a .bam
print STDERR "Converting final putative lgt reads into .bam . . .\n";
`cat $output_prefix.Ff-merged.PE-filtered.sam >> $output_prefix.header`;
`samtools view -bS $output_prefix.header > $output_prefix.Ff-merged.PE-filtered.bam`;


## Pull out UM_UM as potential Bacterial reads.
if ($options{bac}){
	print STDERR "Parsing out UM_UM putative Bacterial reads . . .\n";
    `samtools view -b -F0x600 -f0xC $options{input_bam} > $output_prefix.UM_UM.bam`;
    if ($? != 0 ){die "Fatal Error line 47.\n";}
}

## Cleanup filtering files
if ($options{cleanup}){
    print STDERR "Removing unnecessary intermediate files . . .\n";
    `rm $output_prefix.F4.f8.bam`;
    `rm $output_prefix.F8.f4.bam`;
    `rm $output_prefix.Ff-merged.PE-filtered.sam`;
    `rm $output_prefix.Ff-merged.bam`;
    `rm $output_prefix.header`;
    if ($? != 0 ){die "Fatal Error line 41.\n";}
}



print STDERR "Finished processing $options{input_bam}\n";
