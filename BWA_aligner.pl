#!/usr/bin/perl -I /home/ksieber/scripts/ -I /home/ksieber/perl5/lib/perl5/
use strict;
use warnings;
use File::Basename;
use lib ( '/home/ksieber/scripts/', '/local/projects-t3/HLGT/scripts/lgtseek/lib' );
use run_cmd;
use print_call;
use Cwd;
use bwa;
use setup_input;
use POSIX;

if ( !@ARGV ) { &help; }

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
our $results = GetOptions(
    \%options,    'input|i=s',    'input_list=s',     'output_prefix|p=s', 'output_dir|o=s', 'subdirs=s',
    'ref|r=s',    'ref_list=s',   'threads|t=i',      'disable_SW=i',      'bam_output=i',   'sort_index=i',
    'mpileup=i',  'no_cleanup=i', 'insert_metrics=i', 'mapped_only=i',     'cmd_log=s',      'Qsub|q=i',
    'name=s',     'project=s',    'sub_mem=s',        'help|h',            'help_full|?',    'name_sort_input=i',
    'sort_mem=s', 'sub_mail=s',   'mem=i',            'bwa=s',
) or die "Error: Unrecognized command line option. Please try again.\n";

## Help subroutines (at the end of the sciprt)
if ( $options{help} )      { &help; }
if ( $options{help_full} ) { &help_full; }

## Default values
if ( !$options{input} && !$options{input_list} ) {
    die "Error: Must give input files to map with
 --input or --input_list.\n";
}
if ( !$options{ref} && !$options{ref_list} ) { die "ERROR:  Must have enter a reference file to use.\n"; }
if ( !defined $options{output_dir} ) { $options{output_dir} = getcwd; }
run_cmd("mkdir -p $options{output_dir}");
my @in_suffix_list = ( '.bam', '.fastq.gz', '_\d+.fastq', '.fastq', '.fq' );
my @ref_suffix_list = ( '.fasta', '.fa', '.fna', '.txt' );
$options{sort_index}
    = ( defined $options{mpileup} and $options{mpileup} == 1 )
    ? 1
    : $options{sort_index};    ## Mpileup needs a srt.bam, so turn it on automatically if --mpileup is passed
my $subdirs  = defined $options{subdirs}  ? "$options{subdirs}"  : "0";
my $threads  = defined $options{threads}  ? "$options{threads}"  : "1";
my $project  = defined $options{project}  ? "$options{project}"  : "jdhotopp-lab";
my $sub_mem  = defined $options{sub_mem}  ? "$options{sub_mem}"  : "6G";
my $sort_mem = defined $options{sort_mem} ? "$options{sort_mem}" : "5G";

if ( defined $options{Qsub} and ( $options{'Qsub'} == 1 && defined $options{name_sort_input} and $options{name_sort_input} == 1 ) ) {
    $sub_mem =~ /^(\d+)[K|M|G]$/;    ## remove "G" from sub_mem=#G ie "gigs"
    my $sub_mem_quant = $1;
    $sort_mem =~ /^(\d+)[K|M|G]$/;
    my $sort_mem_quant = $1;

    if ( $sub_mem_quant < $sort_mem_quant * $threads ) {
        $sub_mem = ( ceil( ( $sort_mem_quant * $threads ) * 1.1 ) ) + 1 . "G";
    }
}

print_notebook( \%options );
my @ref_list;

## Setup the reference list
if ( defined $options{ref} ) {
    push( @ref_list, $options{ref} );
}
if ( defined $options{ref_list} ) {
    open( LIST, "<", "$options{ref_list}" ) || die "Error: Can't open reference list because: $!\n";
    while (<LIST>) {
        chomp;
        push( @ref_list, $_ );
    }
    close LIST;
}

## Setup the input list
my $input = setup_input( \%options );

## Qsub or run BWA
foreach my $files (@$input) {
    foreach my $refs (@ref_list) {
        if ( $subdirs == 1 ) {
            my ( $f1, $f2 ) = split( /,/, $files );
            $f1 =~ /\/(\w+)\.\w+$/;
            my $name = defined $options{output_prefix} ? $options{output_prefix} : $1;
            run_cmd("mkdir -p $options{output_dir}");    ## Make sure we have the original output_dir
            $options{output_dir} = "$options{output_dir}/" . "$name/";    ## Add the subdir to the output_dir name
            run_cmd("mkdir -p $options{output_dir}");                     ## Make the subdir if we need to
        }
        if ( defined $options{Qsub} and $options{Qsub} == 1 ) {
            my $cmd = "/home/ksieber/scripts/BWA_aligner.pl";
            if ( $options{input_list} ) {
                $options{input} = $files;
            }    ## If we are in the orignal call, we need to make sure to qsub a single input
            if ( $options{ref_list} ) { $options{ref} = $refs; }
            foreach my $key ( keys %options ) {
                next if ( $options{input_list} && $key =~ /input_list/ );    ## If we are in the orignal call, we don't want to qsub more lists
                next if ( $options{ref_list}   && $key =~ /ref_list/ );
                next if ( $key =~ /subdirs/ );                               ## We already setup the subdirs so we skip it now
                next
                    if ( $key =~ /Qsub/ && !$options{input_list} )
                    ;    ## If we are in the orignal call with input_list, we probably want to qsub each input
                if ( $options{$key} ) { $cmd = $cmd . " --$key=$options{$key}"; }
            }
            my $job_name = defined $options{name} ? "$options{name}" : "BWA_aln";
            Qsub(
                {   cmd      => "$cmd",
                    threads  => "$threads",
                    sub_mem  => "$sub_mem",
                    sub_name => $job_name,
                    sub_mail => $options{sub_mail},
                    wd       => "$options{output_dir}",
                    project  => "$project",
                }
            );
            next;
        }
        elsif ( defined $options{mem} and $options{mem} == 1 ) {
            bwa_mem( $files, $refs, \%options );
        }
        else {
            bwa_aln( $files, $refs, \%options );
        }
    }
}

############################
#### --Help subroutines ####
############################
sub help {
    die "\nHELP: This script will BWA align the input to a reference.
        --input=            Input file to be BWA mapped. Either: in.bam or in_1.fq,in_2.fq
        --ref=              Reference.fna+index
        --output_dir=       Will output to current working directory unless another is specified with this. ie. /ksieber_dir/tmp/
        --help|h            Basic Help Info
        --help_full|?       Full Help Info\n\n";
}

sub help_full {
    die "\nHELP: This script will align the input (fastq/bam) to a reference.
    --input|i=              Input file to be BWA mapped. Either: in.bam or in_1.fq,in_2.fq
    --input_list=           List of input files to be mapped. 1 bam/line. _1,_2 fastq/line (fastqs MUST be comma seperated).
    --name_sort_input=      <0|1> [0] 1=Name sort the input bam
    --ref|r=                Reference.fna+index
    --ref_list=             List of References.
    --output_prefix=        Prefix for each output. Default output = \$prefix_at_\$ref_name.bam
    --output_dir|o=         /directory/for/output/ 
      --subdirs=            <0|1> [0] 1= Make subdirectories for each input file to be mapped. 
    --mem=                  <0|1> [0] 1= Use bwa mem instead of aln.
    --disable_SW=           <0|1> [0] 1= Disable Smith-Waterman for the UM mate. Ideal for quicker LGT mappings IF they are high confidence. 
    --mapped_only=          <0|1> [0] 1= Only keep mates with 1 mapped read.
    --sam_output=           <0|1> [0] 0= .bam output; 1= .sam output
    --sort_index=           <0|1> [0] 1= Sort and index the new.bam into new.srt.bam and new.srt.bai. 
    --mpileup=              <0|1> [0] 1= Calculate pileup coverage on .bam.
    --no_cleanup=           <0|1> [0] 0= Removes .sai files and unsorted.bam with --sort_index. 1=No deleting intermediate data. 
    --insert_metrics=       <0|1> [0] 1= Use Picard to calculate insert size metrics.
    --Qsub|q=               <0|1> [0] 1= qsub the mapping to SGE grid.
      --threads|t=          < # >   [1] Set the number of cpu threads to use for bwa aln steps. USE CAREFULLY.
      --project=            [jdhotopp-lab].
      --sub_mem=            [6G] Memory free for qsub.
      --sub_mail=           [0] 1= email user\@som.umaryland.edu when job is complete & with stats. Can also specify --sub_mail=specific\@email.foo
      --sub_name=           Name qsub submission.
      --wd=                 [--output_dir]
    --cmd_log=              <0|1> [0] 1= Log all commands run in each output_dir/output_prefix.cmd_log  
    --help\n\n";
}

__END__
