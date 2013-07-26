#!/usr/local/bin/perl -w
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
                          'ref_list=s',
                          'fastq_1=s',
                          'fastq_2=s',
                          'bam_input=s',
                          'output_dir=s',
                          't=s',
                          'nodes=s',
                          'bam_output',
                          'sort_index_bams',
                          'mpileup',
                          'cleanup',
                          'insert_metrics',
                          'help|h',
                         );

if ($options{help}){die "\nHELP: This script will take a list of references and make shell scripts align the input (fastq/bam) to ref.list 
\t--ref=               REQUIRED.   Reference.fna+index
\t--fastq_1=	       ~Required.  fastq_1 file
\t--fastq_2=	       ~Required.  fastq_2 file
\t--bam_input=         ~Required.  data.bam
\t--fastq_list=        List of fastq files to align.  1 pair of _1&_2 per line. Each line is a different fastq pair.        
\t--output_prefix=     REQUIRED.   Prefix for each output.  Ie. SRA_LGT_NC001234
\t--output_dir=        REQUIRED.   Directory for output. ie /foo/bar/tmp/
\t--t=                 # [1] Will set the number of cpu threads to use for bwa aln steps. USE CAREFULLY.   
\t--nodes=             # [150] This will set the number of shell scripts to make. 
\t--bam_output         Convert the new.sam into a new.bam
\t--sort_index_bams    Sort and index the new.bam into new.srt.bam and new.srt.bai.  If --bam_output wasn't used it will generate the .bam anways. 
\t--mpileup            Calculate pileup coverage on .bam.  Must be used with --bam_output
\t--cleanup            WARNING: This option delete new.1.sai; new.2.sai; and w/ --bam_output it will rm new.sam; new.bam Leaving only new.srt.bam/bai
\t--insert_metrics     Use Picard to calculate insert size metrics. Must be used with --bam_output.
\t--help
";}

## Make sure all the proper inputs are given.
if (!$options{fastq_1} && !$options{fastq_2} && !$options{bam_input}){die "ERROR: Must give an input to align.  Use either --fastq_1/_2, or bam_inputn.\n";}
if (!$options{ref_list}){die "ERROR:  Must have enter a reference file to use.\n";}
if (!$options{output_prefix} || !$options{output_dir}){die "ERROR: Must use --output_prefix= and --output_dir=.\n";}
if ($options{insert_metrics} && !$options{bam_output}){die "ERROR: Must use --bam_output with --insert_metrics.\n";}
my $threads = $options{t} ? "$options{t}" : "1";
my $nodes = $options{nodes} ? "$options{nodes}" : "150";

open (LIST, "<", $options{ref_list}) || die "ERROR: Can't open $options{ref_list} because: $!\n";
my @ref_list;
while(<LIST>){
    chomp;
    push(@ref_list, $_);
}
close LIST || die "ERROR: Can't close $options{ref_list} because: $!\n";

my $ref_list_length = @ref_list;
my $shell_script_size = ($ref_list_length/$nodes)+1;
my $ref_counter = 0;
for (my $i=0; $i<$nodes; $i++){
    open(SHELL, ">", "$options{output_dir}/Refloop_$i\.sh");
    for (my $f=0; $f<$shell_script_size; $f++){
        print SHELL "perl /local/projects/HLGT/ksieber_dir/perl/BWA_aligner.pl --ref=$ref_list[$ref_counter] --output_prefix=$options{output_prefix} --output_dir=$options{output_dir} --t=$threads";
        if ($options{fastq_1}){print SHELL " --fastq_1=$options{fastq_1} --fastq_2=$options{fastq_2}";}
        if ($options{bam_input}){print SHELL " --bam_input=$options{bam_input}";}
        if ($options{bam_output}){print SHELL " --bam_output";}
        if ($options{sort_index_bams}){print SHELL " --sort_index_bams";}
        if ($options{mpileup}){print SHELL " --mpileup";}
        if ($options{cleanup}){print SHELL " --cleanup";}
        print SHELL "\n";
        $ref_counter++;
    }
    close SHELL || die "ERROR: Can't close $options{output_dir}/Refloop_$i\.sh because: $!\n";
}
