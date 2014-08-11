#!/usr/bin/perl -I /home/ksieber/scripts/ -I /home/ksieber/perl5/lib/perl5/
use strict;
use File::Basename;
use lib ( '/home/ksieber/scripts/', '/local/projects-t3/HLGT/scripts/lgtseek/lib' );
use run_cmd;
use setup_input;
use POSIX;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
our $results = GetOptions(
    \%options,     'input|i=s',    'input_list=s', 'output_prefix=s',     'output_dir|o=s', 'subdirs=s',    'ref|r=s',          'ref_list=s',
    'threads|t=s', 'disable_SW=s', 'bam_output=s', 'sort_index_output=s', 'mpileup=s',      'no_cleanup=s', 'insert_metrics=s', 'mapped_only=s',
    'cmd_log=s',   'Qsub|Q=i',     'name=s',       'project=s',           'sub_mem=s',      'help|h',       'help_full',        'name_sort_input=i',
    'sort_mem=s',  'sub_mail=s',
) or die "Error: Unrecognized command line option. Please try again.\n";

## Help subroutines (at the end of the sciprt)
if ( $options{help} )      { &help; }
if ( $options{help_full} ) { &help_full; }

## Default values
if ( !$options{input} && !$options{input_list} ) {
    die "Error: Must give input files to map with --input or --input_list.\n";
}
if ( !$options{output_dir} ) { die "Error: Must use --output_dir=/path/to/output/\n"; }
run_cmd("mkdir -p $options{output_dir}");
if ( !$options{ref} && !$options{ref_list} ) { die "ERROR:  Must have enter a reference file to use.\n"; }
my @in_suffix_list = ( '.bam', '.fastq.gz', '_\d+.fastq', '.fastq', '.fq' );
my @ref_suffix_list = ( '.fasta', '.fa', '.fna', '.txt' );
my $threads = defined $options{t} ? "$options{t}" : "1";    ## Default # of threads = 1
$options{sort_index_output} = ( $options{mpileup} == 1 ) ? "1" : "$options{sort_index_output}";    ## Mpileup needs a srt.bam, so turn it on automatically if --mpileup is passed
my $subdirs  = defined $options{subdirs}  ? "$options{subdirs}"  : "0";
my $threads  = defined $options{threads}  ? "$options{threads}"  : "1";
my $project  = defined $options{project}  ? "$options{project}"  : "jdhotopp-lab";
my $sub_mem  = defined $options{sub_mem}  ? "$options{sub_mem}"  : "6G";
my $sort_mem = defined $options{sort_mem} ? "$options{sort_mem}" : "5G";

if ( $options{'Qsub'} == 1 && $options{name_sort_input} == 1 ) {
    $sub_mem =~ /^(\d+)[K|M|G]$/;                                                                  ## remove "G" from sub_mem=#G ie "gigs"
    my $sub_mem_quant = $1;
    $sort_mem =~ /^(\d+)[K|M|G]$/;
    my $sort_mem_quant = $1;

    if ( $sub_mem_quant < $sort_mem_quant * $threads ) {
        $sub_mem = ( ceil( ( $sort_mem_quant * $threads ) * 1.1 ) ) + 1 . "G";
    }
}

my @ref_list;

## Setup the reference list
if ( $options{ref} ) {
    push( @ref_list, $options{ref} );
}
if ( $options{ref_list} ) {
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
        if ( $options{Qsub} == 1 ) {
            my $cmd = "/home/ksieber/scripts/BWA_aligner.pl";
            if ( $options{input_list} ) {
                $options{input} = $files;
            }                                                             ## If we are in the orignal call, we need to make sure to qsub a single input
            if ( $options{ref_list} ) { $options{ref} = $refs; }
            foreach my $key ( keys %options ) {
                next
                    if ( $options{input_list} && $key =~ /input_list/ );    ## If we are in the orignal call, we don't want to qsub more lists
                next if ( $options{ref_list} && $key =~ /ref_list/ );
                next if ( $key =~ /subdirs/ );                              ## We already setup the subdirs so we skip it now
                next
                    if ( $key =~ /Qsub/ && !$options{input_list} );         ## If we are in the orignal call with input_list, we probably want to qsub each input
                if ( $options{$key} ) { $cmd = $cmd . " --$key=$options{$key}" }
                ;                                                           ## Build the command for all other options passed in @ original call
            }
            my $job_name = defined $options{name} ? "$options{name}" : "BWA_aln";
            Qsub(
                {   cmd      => "$cmd",
                    threads  => "$threads",
                    mem      => "$sub_mem",
                    wd       => "$options{output_dir}",
                    sub_name => $job_name,
                    project  => "$project",
                }
            );
            next;
        }
        else {
            bwa_align( $files, $refs );
        }
    }
}

## BWA alignment
sub bwa_align {
    my ( $files, $ref ) = @_;
    print "@_\n";
    my ( $ref_name, $ref_path, $ref_suf ) = fileparse( $ref, @ref_suffix_list );    ## Grab the reference name to use for the naming the output
    my $file1;                                                                      ## Global input file name
    my $file2;                                                                      ## Global input file name2
    my $bam;                                                                        ## 0=fastq, 1=bam
    if ( $files =~ /.fq$/ || $files =~ /.fastq$/ || $files =~ /.fastq.gz$/ ) {
        if ( $options{name_sort_input} == 1 ) { die "Error: This script can't name-sort fastq files. Please adjust & try again.\n"; }
        ( $file1, $file2 ) = split( /,/, $files );
        if ( $file2 !~ /\w+/ ) { $file2 = $file1; }
        $bam = 0;
    }
    elsif ( $files =~ /.bam$/ ) {
        $file1 = $files;
        $file2 = $files;
        $bam   = 1;
    }
    else {
        die "Could not resolve the input type. Make sure it is either a .bam,.fq,.fastq\n";
    }

    ## Setup log
    my $log;
    if ( $options{cmd_log} == 1 ) {
        $log = "$options{output_dir}/log.txt";
    }

    ## Setup output prefix (path/file-name).bam
    my ( $input, $path, $suf ) = fileparse( $file1, @in_suffix_list );
    my $out = $options{output_prefix} ? "$options{output_prefix}" : "$input\_at_$ref_name";
    my $dir = $options{output_dir};
    if ( $dir =~ /\/\/$/ ) { $dir =~ s/\/$//g; }
    my $output_prefix = "$dir\/$out";

    ## Name sort input
    if ( $files =~ /.bam$/ && $options{name_sort_input} == 1 ) {
        run_cmd( "samtools sort -n -@ $threads -m $sort_mem $file1 $output_prefix\_name-sorted", $log );
        $file1 = "$output_prefix\_name-sorted.bam";
        $file2 = "$output_prefix\_name-sorted.bam";
    }

    ## Run BWA ALN
    if ( $bam == 1 ) {
        run_cmd( "bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref -b1 $file1 > $output_prefix\.1.sai 2>>$output_prefix\_bwa_stderr.log", $log );
        run_cmd( "bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref -b2 $file2 > $output_prefix\.2.sai 2>>$output_prefix\_bwa_stderr.log", $log );
    }
    else {
        run_cmd( "bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref $file1 > $output_prefix\.1.sai 2>>$output_prefix\_bwa_stderr.log", $log );
        run_cmd( "bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref $file2 > $output_prefix\.2.sai 2>>$output_prefix\_bwa_stderr.log", $log );
    }

    ## BWA SAMPE
    if ( $options{mapped_only} == 1 ) {
        run_cmd( "bwa sampe $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 | samtools view -F0x4 -bhS - > $output_prefix\.bam", $log );
    }
    elsif ( $options{sam_output} == 1 ) {
        if ( $options{disable_SW} == 1 ) {
            run_cmd( "bwa sampe -s $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 > $output_prefix\.sam 2>>$output_prefix\_bwa_stderr.log", $log );
        }
        else {
            run_cmd( "bwa sampe $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 > $output_prefix\.sam 2>>$output_prefix\_bwa_stderr.log", $log );
        }
    }
    else {
        if ( $options{disable_SW} == 1 ) {
            run_cmd( "bwa sampe -s $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 2>>$output_prefix\_bwa_stderr.log | samtools view -bhS - > $output_prefix\.bam ", $log );
        }
        else {
            run_cmd( "bwa sampe $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 2>>$output_prefix\_bwa_stderr.log | samtools view -bhS - > $output_prefix\.bam ", $log );
        }
    }

    ## Sort and index bams
    if ( $options{sort_index_output} == 1 ) {
        run_cmd( "samtools sort $output_prefix\.bam $output_prefix\.srt",          $log );
        run_cmd( "samtools index $output_prefix\.srt.bam $output_prefix\.srt.bai", $log );
    }

    ## Samtools Mpileup
    if ( $options{mpileup} == 1 ) {
        run_cmd( "samtools mpileup -Af $options{ref} $output_prefix\.srt.bam > $output_prefix\.COVERAGE.txt", $log );
    }

    ## Picard insert metrics
    if ( $options{insert_metrics} == 1 ) {
        run_cmd(
            "java -Xmx3g -jar /home/jdhotopp/bin/Picard/picard-tools-1.48/CollectInsertSizeMetrics.jar AS=true I=$output_prefix.srt.bam O=$output_prefix\_std_insert.metrics H=$output_prefix\_std_insert.histogram M=0 VALIDATION_STRINGENCY=SILENT",
            $log
        );
        run_cmd(
            "java -Xmx3g -jar /home/jdhotopp/bin/Picard/picard-tools-1.48/CollectInsertSizeMetrics.jar AS=true I=$output_prefix.srt.bam O=$output_prefix\_lrg_insert.metrics H=$output_prefix\_lrg_insert.histogram M=0 VALIDATION_STRINGENCY=SILENT DEVIATIONS=1000000000000000000",
            $log
        );
    }

    ## Cleanup intermediate files (.sai)
    unless ( $options{no_cleanup} == 1 ) {
        run_cmd( "rm $output_prefix\.1.sai", $log );
        run_cmd( "rm $output_prefix\.2.sai", $log );
        if ( $options{sort_index_output} == 1 ) {
            run_cmd( "rm $output_prefix\.bam", $log );
        }
        run_cmd( "rm $output_prefix\_bwa_stderr.log", $log );
    }
    print STDERR "====== Completed BWA mapping: $file1 against: $ref output: $output_prefix ======\n";
}

#### --Help subroutines ####
############################
sub help {
    die "\nHELP: This script will BWA align the input to a reference.
        --input=            Input file to be BWA mapped. Either: in.bam or in_1.fq,in_2.fq
        --ref=              Reference.fna+index
        --output_dir=           Will output to current working directory unless another is specified with this. ie. /ksieber_dir/tmp/
        --help              Basic Help Info
        --help_full             Full Help Info\n";
}

sub help_full {
    die "\nHELP: This script will align the input (fastq/bam) to a reference.
    --input|i=              Input file to be BWA mapped. Either: in.bam or in_1.fq,in_2.fq
    --input_list=           List of input files to be mapped. 1 bam/line. _1,_2 fastq/line (fastqs MUST be comma seperated).
    --ref|r=                Reference.fna+index
    --ref_list=             List of References.
    --output_prefix=        Prefix for each output.  Ie. (SRA_LGT)_at_\$ref_name
    --output_dir|o=         Will output to current working directory unless another is specified with this. ie. /ksieber_dir/tmp/
    --subdirs=              <0|1> [0] 1= Make subdirectories for each input file to be mapped. 
    --disable_SW=           <0|1> [0] 1= Disable Smith-Waterman for the UM mate. Ideal for quicker LGT mappings IF they are high confidence. 
    --mapped_only=          <0|1> [0] 1= Only keep mates with 1 mapped read.
    --sam_output=           <0|1> [0] 0= .bam output; 1= .sam output
    --sort_index_output=    <0|1> [0] 1= Sort and index the new.bam into new.srt.bam and new.srt.bai. 
    --mpileup=              <0|1> [0] 1= Calculate pileup coverage on .bam.
    --no_cleanup=           <0|1> [0] 0= Removes .sai files and unsorted.bam with --sort_index_output. 1=No deleting intermediate data. 
    --insert_metrics=       <0|1> [0] 1= Use Picard to calculate insert size metrics.
    --Qsub|Q=               <0|1> [0] 1= qsub the mapping to SGE grid.
      --threads|t=          < # >   [1] Set the number of cpu threads to use for bwa aln steps. USE CAREFULLY.
      --project=            [jdhotopp-lab].
      --sub_mem=            [6G] Memory free for qsub.
      --sub_mail=           [0] 1= email user\@som.umaryland.edu when job is complete & with stats. Can also specify --sub_mail=specific\@email.foo
      --sub_name=           Name qsub submission.
      --wd=                 [--output_dir]
    --cmd_log=              <0|1> [0] 1= Log all commands run in each output_dir/output_prefix.cmd_log  
    --help\n";
}

__END__
