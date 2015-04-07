package bwa;
use warnings;
use strict;
use run_cmd;
use mk_dir;
use File::Basename;
use Carp;
$Carp::MaxArgLen = 0;    ## Report full length error
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw( bwa_aln bwa_mem );

=head1

Title   : bwa_mem
Usage   : my $mapped_bam = bwa_mem($input,$reference,{other_options => possible})
Function: BWA mem an input file against a reference.
Args    : 
    input=          *Mandatory* Input file to be BWA mapped. Either: in.bam or in_1.fq. If in_1.fq and no in_2.fq, will auto look for in_2.fq. can also use =in_1.fq,in_2.fq
    ref=            *Mandatory* Reference.fna+index
    {
        output          =>  /full/path/and/name.bam for the output (overrides output_dir & output_prefix)
        output_dir      =>  Will output to current working directory unless another is specified with this. ie. /ksieber_dir/tmp/
        output_prefix   =>  Prefix for each output.  ex out= $output_dir/($prefix)_mem_($ref_name)
        sort_index      =>  <0|1> [0] 1= Sort and index the new.bam into new.srt.bam and new.srt.bai. 
        sort_mem        =>  [1G] RAM / cpu to sort by.
        threads         =>  [1]  Number of cpu to use. 
        interleaved     =>  <0|1> [0] 1= first query file consists of interleaved paired-end sequences.
        other           =>  A string with other bwa mem options. ex: other => "-j chr1",
        Picard_jar      =>  [/home/ksieber/lib/picard/dist/picard.jar] Path to: picard.jar. Only needed if input=bam or insert_metrics=1.
        bwa             =>  [/home/ksieber/lib/bwa/bwa] Path to bwa executable. 
        mpileup         =>  <0|1> [0] 1= Calculate pileup coverage on .bam.
        insert_metrics  =>  <0|1> [0] 1= Use Picard to calculate insert size metrics.
        no_cleanup      =>  <0|1> [0] 0= Removes created fq files and unsorted.bam if --sort_index=1. 1=No deleting intermediate data. 
    }

Returns : File path to a newly mapped Bam.

=cut

sub bwa_mem {
    my ( $input, $ref, $options ) = @_;

    my @in_suffix_list = ( '.bam', '_\d+.fastq', '.fastq', '.fq', '.fastq.gz', );
    my @ref_suffix_list = ( '.fasta', '.fa', '.fna', '.txt' );
    my $threads    = defined $options->{threads}    ? $options->{threads}    : "1";
    my $sort_mem   = defined $options->{sort_mem}   ? $options->{sort_mem}   : "1G";
    my $other      = defined $options->{other}      ? "$options->{other} "   : "";
    my $Picard_jar = defined $options->{Picard_jar} ? $options->{Picard_jar} : "/home/ksieber/lib/picard/dist/picard.jar";
    my $bwa        = defined $options->{bwa}        ? $options->{bwa}        : "/home/ksieber/lib/bwa/bwa";
    my $interleaved = ( defined $options->{interleaved} and $options->{interleaved} == 1 ) ? "-p " : "";
    my $fastq1;    ## Global input file name
    my $fastq2;    ## Global input file name2

    ## Setup output prefix (path/file-name).bam
    my ( $ref_name, $ref_path, $ref_suf ) = fileparse( $ref,   @ref_suffix_list );    ## Grab the reference name to use for the naming the output
    my ( $in_fn,    $in_path,  $in_suf )  = fileparse( $input, @in_suffix_list );
    my $out_prefix = $options->{output_prefix} ? $options->{output_prefix} : "$in_fn\_mem_$ref_name";
    my $out_dir    = $options->{output_dir}    ? $options->{output_dir}    : $in_path;
    if ( defined $options->{output} ) { ( $out_prefix, $out_dir, $ref_suf ) = fileparse( $options->{output}, '.bam' ); }
    if ( $out_dir =~ /\/\/$/ ) { $out_dir =~ s/\/$//g; }
    mk_dir($out_dir);
    my $output_dir_pref = "$out_dir/$out_prefix";

    if ( $in_suf eq ".bam" ) {
        if ( !-e $input ) { confess "Error: input bam(?) does not exist: $input"; }
        run_cmd("java -Xmx4G -jar $Picard_jar SamToFastq VALIDATION_STRINGENCY=SILENT I=$input  F=$output_dir_pref\_1.fastq F2=$output_dir_pref\_2.fastq");
        $fastq1 = "$output_dir_pref\_1.fastq";
        $fastq2 = "$output_dir_pref\_2.fastq";
    }
    else {
        ( $fastq1, $fastq2 ) = split( /,/, $input );
        if ( !-e $fastq1 ) { confess "Error: fastq1 does not exist: $fastq1\n"; }
        if ( !$fastq2 ) {
            my ( $fq_fn, $fq_path, $fq_suf ) = fileparse( $fastq1, ( "_1.fastq", "_1.fq" ) );
            if    ( -e "$fq_path/$fq_fn\_2.fastq" ) { $fastq2 = "$fq_path/$fq_fn\_2.fastq"; }
            elsif ( -e "$fq_path/$fq_fn\_2.fq" )    { $fastq2 = "$fq_path/$fq_fn\_2.fq"; }
            else                                    { $fastq2 = ""; }
        }

    }

    ## Setup log
    my $log;
    if ( $options->{cmd_log} ) {
        $log = "$out_dir/log.txt";
    }

    run_cmd( "$bwa mem $interleaved$other-t $threads $ref $fastq1 $fastq2 | samtools view -S - -bo $output_dir_pref\.bam", $log );

    ## Sort and index bams
    if ( ( defined $options->{sort_index} and $options->{sort_index} == 1 ) or ( defined $options->{insert_metrics} and $options->{insert_metrics} == 1 ) ) {
        run_cmd( "samtools sort -@ $threads -m $sort_mem $output_dir_pref\.bam $output_dir_pref\.srt 2>>$output_dir_pref\_bwa_stderr.log", $log );
        run_cmd( "samtools index $output_dir_pref\.srt.bam 2>>$output_dir_pref\_bwa_stderr.log",                                           $log );
    }

    ## Samtools Mpileup
    if ( defined $options->{mpileup} and $options->{mpileup} == 1 ) {
        run_cmd( "samtools mpileup -Af $ref $output_dir_pref\.srt.bam > $output_dir_pref\.COVERAGE.txt 2>>$output_dir_pref\_bwa_stderr.log", $log );
    }

    ## Picard insert metrics
    if ( defined $options->{insert_metrics} and $options->{insert_metrics} == 1 ) {
        run_cmd(
            "java -Xmx4g -jar $Picard_jar CollectInsertSizeMetrics AS=true I=$output_dir_pref\.srt.bam O=$output_dir_pref\_std_insert.metrics H=$output_dir_pref\_std_insert.histogram M=0 VALIDATION_STRINGENCY=LENIENT",
            $log
        );
        run_cmd(
            "java -Xmx4g -jar $Picard_jar CollectInsertSizeMetrics AS=true I=$output_dir_pref\.srt.bam O=$output_dir_pref\_lrg_insert.metrics H=$output_dir_pref\_lrg_insert.histogram M=0 VALIDATION_STRINGENCY=LENIENT DEVIATIONS=1000000000000000000",
            $log
        );
    }

    ## Cleanup intermediate files (.sai)
    unless ( defined $options->{no_cleanup} and $options->{no_cleanup} == 1 ) {
        if ( $in_suf eq ".bam" ) {
            run_cmd( "rm $fastq1", $log );
            run_cmd( "rm $fastq2", $log );
        }
        if ( ( defined $options->{sort_index} and $options->{sort_index} == 1 ) or ( defined $options->{insert_metrics} and $options->{insert_metrics} == 1 ) ) {
            run_cmd( "rm $output_dir_pref\.bam", $log );
        }
        if ( -e "$output_dir_pref\_bwa_stderr.log" ) { run_cmd( "rm $output_dir_pref\_bwa_stderr.log", $log ); }
    }
    my $retval = defined $options->{sort_index} ? "$output_dir_pref\.srt.bam" : "$output_dir_pref\.bam";

    # print STDERR "====== Completed BWA mapping: $file1 against: $ref output: $output_prefix\.bam ======\n";
    return $retval;
}

=head1

Title   : bwa_aln
Usage   : my $mapped_bam = bwa_aling($input,$reference,{other_options => possible})
Function: BWA align an input file against a reference.
Args    : 
    input=          *Mandatory* Input file to be BWA mapped. Either: in.bam or in_1.fq,in_2.fq
    ref=            *Mandatory* Reference.fna+index
    {
        output_prefix   =>  Prefix for each output.  Ie. (SRA_LGT)_at_\$ref_name
        output_dir      =>  Will output to current working directory unless another is specified with this. ie. /ksieber_dir/tmp/
        sort_index      =>  <0|1> [0] 1= Sort and index the new.bam into new.srt.bam and new.srt.bai. 
        M_only          =>  <0|1> [0] 1= Only keep mates with M_*
        MM_only         =>  <0|1> [0] 1= Only keep mates with M_M
        mpileup         =>  <0|1> [0] 1= Calculate pileup coverage on .bam.
        insert_metrics  =>  <0|1> [0] 1= Use Picard to calculate insert size metrics.
        disable_SW      =>  <0|1> [0] 1= Disable Smith-Waterman for the UM mate. Ideal for quicker LGT mappings IF they are high confidence. 
        sam_output      =>  <0|1> [0] 0= .bam output; 1= .sam output
        no_cleanup      =>  <0|1> [0] 0= Removes .sai files and unsorted.bam with --sort_index. 1=No deleting intermediate data. 
    }

Returns : File path to a newly mapped Bam.

=cut

sub bwa_aln {
    my ( $files, $ref, $options ) = @_;

    my @in_suffix_list = ( '.bam', '.fastq.gz', '_\d+.fastq', '.fastq', '.fq' );
    my @ref_suffix_list = ( '.fasta', '.fa', '.fna', '.txt' );
    my $threads    = defined $options->{threads}    ? $options->{threads}    : "1";
    my $sort_mem   = defined $options->{sort_mem}   ? $options->{sort_mem}   : "1G";
    my $Picard_jar = defined $options->{Picard_jar} ? $options->{Picard_jar} : "/home/ksieber/lib/picard/dist/picard.jar";
    my $bwa        = defined $options->{bwa}        ? $options->{bwa}        : "bwa";

    ## Setup output prefix (path/file-name).bam
    my ( $ref_name, $ref_path, $ref_suf ) = fileparse( $ref, @ref_suffix_list );    ## Grab the reference name to use for the naming the output
    my $file1;                                                                      ## Global input file name
    my $file2;                                                                      ## Global input file name2
    my $bam;                                                                        ## 0=fastq, 1=bam
    if ( $files =~ /.fq$/ || $files =~ /.fastq$/ || $files =~ /.fastq.gz$/ ) {
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
    if ( $options->{cmd_log} ) {
        $log = "$options->{output_dir}/log.txt";
    }

    ## Setup output prefix (path/file-name).bam
    my ( $input, $path, $suf ) = fileparse( $file1, @in_suffix_list );
    my $out = $options->{output_prefix} ? "$options->{output_prefix}" : "$input\_at_$ref_name";
    my $dir = $options->{output_dir};
    if ( $dir =~ /\/\/$/ ) { $dir =~ s/\/$//g; }
    my $output_prefix = "$dir\/$out";

    if ( $files =~ /\.bam$/ and ( ( defined $options->{name_sort_input} and $options->{name_sort_input} == 1 ) or `samtools view -H $file1 | head -n 1` !~ /queryname/ ) ) {
        $options->{name_sort_input} = 1;
        run_cmd( "samtools sort -n -@ $threads -m $sort_mem $file1 $output_prefix\_name-sorted", $log );
        $file1 = "$output_prefix\_name-sorted.bam";
        $file2 = "$output_prefix\_name-sorted.bam";
    }

    ## Run BWA ALN
    if ( $bam == 1 ) {
        run_cmd( "$bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref -b1 $file1 > $output_prefix\.1.sai 2>>$output_prefix\_bwa_stderr.log", $log );
        run_cmd( "$bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref -b2 $file2 > $output_prefix\.2.sai 2>>$output_prefix\_bwa_stderr.log", $log );
    }
    else {
        run_cmd( "$bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref $file1 > $output_prefix\.1.sai 2>>$output_prefix\_bwa_stderr.log", $log );
        run_cmd( "$bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref $file2 > $output_prefix\.2.sai 2>>$output_prefix\_bwa_stderr.log", $log );
    }

    if ( defined $options->{M_only} and $options->{M_only} ) {
        run_cmd( "$bwa sampe $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 | samtools view -bhS -F0x4 - > $output_prefix\.bam 2>>$output_prefix\_bwa_stderr.log", $log );
    }
    elsif ( defined $options->{MM_only} and $options->{MM_only} ) {
        run_cmd( "$bwa sampe $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 | samtools view -bhS -F0xC - > $output_prefix\.bam 2>>$output_prefix\_bwa_stderr.log", $log );
    }
    elsif ( defined $options->{sam_output} and $options->{sam_output} == 1 ) {
        if ( defined $options->{disable_SW} and $options->{disable_SW} == 1 ) {
            run_cmd( "$bwa sampe -s $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 > $output_prefix\.sam 2>>$output_prefix\_bwa_stderr.log", $log );
        }
        else {
            run_cmd( "$bwa sampe $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 > $output_prefix\.sam 2>>$output_prefix\_bwa_stderr.log", $log );
        }
    }
    else {
        if ( defined $options->{disable_SW} and $options->{disable_SW} == 1 ) {
            run_cmd(
                "$bwa sampe -s $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 2>>$output_prefix\_bwa_stderr.log | samtools view -bhS - > $output_prefix\.bam 2>>$output_prefix\_bwa_stderr.log",
                $log
            );
        }
        else {
            run_cmd(
                "$bwa sampe $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 2>>$output_prefix\_bwa_stderr.log | samtools view -bhS - > $output_prefix\.bam 2>>$output_prefix\_bwa_stderr.log",
                $log
            );    ## 2>>$output_prefix\_bwa_stderr.log
        }
    }

    ## Sort and index bams
    if (   ( defined $options->{sort_index} and $options->{sort_index} == 1 )
        or ( defined $options->{mpileup} and $options->{mpileup} == 1 )
        or ( defined $options->{insert_metrics} and $options->{insert_metrics} == 1 ) )
    {
        run_cmd( "samtools sort $output_prefix\.bam $output_prefix\.srt 2>>$output_prefix\_bwa_stderr.log", $log );
        run_cmd( "samtools index $output_prefix\.srt.bam 2>>$output_prefix\_bwa_stderr.log",                $log );
    }

    ## Samtools Mpileup
    if ( defined $options->{mpileup} and $options->{mpileup} == 1 ) {
        if ( $ref !~ /\w+\.f\w*a$/ ) {    # Check if the ref. is an actual .fasta file. If not, calc mpileup w/o ref.
            run_cmd( "samtools mpileup -A $output_prefix\.srt.bam > $output_prefix\.COVERAGE.txt 2>>$output_prefix\_bwa_stderr.log", $log );
        }
        else {
            run_cmd( "samtools mpileup -Af $ref $output_prefix\.srt.bam > $output_prefix\.COVERAGE.txt 2>>$output_prefix\_bwa_stderr.log", $log );
        }
    }

    ## Picard insert metrics
    if ( defined $options->{insert_metrics} and $options->{insert_metrics} == 1 ) {
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
    if ( !$options->{no_cleanup} or $options->{no_cleanup} != 1 ) {
        run_cmd( "rm $output_prefix\.1.sai",          $log );
        run_cmd( "rm $output_prefix\.2.sai",          $log );
        run_cmd( "rm $output_prefix\_bwa_stderr.log", $log );
        if ( $file1 =~ /\.bam$/ and ( defined $options->{name_sort_input} and $options->{name_sort_input} == 1 ) ) {
            run_cmd( "rm $output_prefix\_name-sorted\.bam", $log );
        }
        if ( defined $options->{sort_index_output} and $options->{sort_index_output} == 1 ) {
            run_cmd( "rm $output_prefix\.bam", $log );
        }
    }
    my $retval = defined $options->{sort_index} ? "$output_prefix\.srt.bam" : "$output_prefix\.bam";

    # print STDERR "====== Completed BWA mapping: $file1 against: $ref output: $output_prefix\.bam ======\n";
    return $retval;
}

1
