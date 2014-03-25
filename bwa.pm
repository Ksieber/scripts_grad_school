package bwa;
use warnings;
use strict;
use run_cmd;
use File::Basename;
use Carp;
$Carp::MaxArgLen = 0; ## Report full length error 
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( bwa_align );

=head1

Title   : bwa_align
Usage   : my $mapped_bam = bwa_aling($input,$reference,{other_options => possible})
Function: BWA align an input file against a reference.
Args    : 
	input=			*Mandatory* Input file to be BWA mapped. Either: in.bam or in_1.fq,in_2.fq
	ref=			*Mandatory* Reference.fna+index
	{
		output_prefix	=>	Prefix for each output.  Ie. (SRA_LGT)_at_\$ref_name
		output_dir		=>	Will output to current working directory unless another is specified with this. ie. /ksieber_dir/tmp/
		sort_index 		=>	<0|1> [0] 1= Sort and index the new.bam into new.srt.bam and new.srt.bai. 
		M_only 			=>	<0|1> [0] 1= Only keep mates with M_*
		MM_only			=>	<0|1> [0] 1= Only keep mates with M_M
		mpileup 		=>	<0|1> [0] 1= Calculate pileup coverage on .bam.
		insert_metrics 	=> 	<0|1> [0] 1= Use Picard to calculate insert size metrics.
		disable_SW		=>	<0|1> [0] 1= Disable Smith-Waterman for the UM mate. Ideal for quicker LGT mappings IF they are high confidence. 
		sam_output		=>	<0|1> [0] 0= .bam output; 1= .sam output
		no_cleanup 		=>	<0|1> [0] 0= Removes .sai files and unsorted.bam with --sort_index. 1=No deleting intermediate data. 
	}

Returns : File path to a newly mapped Bam.

=cut

sub bwa_align {
	my ($files,$ref,$options) = @_;
	
	my @in_suffix_list=('.bam','.fastq.gz','_\d+.fastq','.fastq','.fq');  
	my @ref_suffix_list=('.fasta','.fa','.fna','.txt');
	my $threads = defined $options->{t} ? "$options->{t}" : "1";

	my ($ref_name,$ref_path,$ref_suf)=fileparse($ref,@ref_suffix_list);     ## Grab the reference name to use for the naming the output
	my $file1;																## Global input file name
	my $file2;																## Global input file name2
	my $bam;     															## 0=fastq, 1=bam
	if($files=~/.fq$/ || $files=~/.fastq$/ || $files=~/.fastq.gz$/){
		($file1,$file2)=split(/,/,$files);
		if($file2!~/\w+/){$file2=$file1;}
		$bam=0;
		} elsif ($files=~/.bam$/){
			$file1=$files;
			$file2=$files;
			$bam=1;
			} else {
				die "Could not resolve the input type. Make sure it is either a .bam,.fq,.fastq\n";
			}

	## Setup log
	my $log;
	if($options->{cmd_log}){
		$log="$options->{output_dir}/log.txt";
	}
	
	## Setup output prefix (path/file-name).bam
	my ($input,$path,$suf)=fileparse($file1,@in_suffix_list);
	my $out = $options->{output_prefix} ? "$options->{output_prefix}" : "$input\_at_$ref_name";
	my $dir = $options->{output_dir};
	if ($dir=~/\/\/$/){$dir =~s/\/$//g;}
	my $output_prefix= "$dir\/$out";


	## Run BWA ALN
	if ($bam==1) {
		run_cmd("bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref -b1 $file1 > $output_prefix\.1.sai 2>>$output_prefix\_bwa_stderr.log",$log);
		run_cmd("bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref -b2 $file2 > $output_prefix\.2.sai 2>>$output_prefix\_bwa_stderr.log",$log);
	} else {
		run_cmd("bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref $file1 > $output_prefix\.1.sai 2>>$output_prefix\_bwa_stderr.log",$log);
		run_cmd("bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref $file2 > $output_prefix\.2.sai 2>>$output_prefix\_bwa_stderr.log",$log);
	}

	## BWA SAMPE
	if ($options->{M_only}){
		run_cmd("bwa sampe $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 | samtools view -F0x4 -bhS - > $output_prefix\.bam",$log);
	} elsif ($options->{MM_only}){
		run_cmd("bwa sampe $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 | samtools view -F0xC -bhS - > $output_prefix\.bam",$log);
	} elsif ($options->{sam_output}) {
		if ($options->{disable_SW}) {
			run_cmd("bwa sampe -s $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 > $output_prefix\.sam 2>>$output_prefix\_bwa_stderr.log",$log);
		} else {
			run_cmd("bwa sampe $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 > $output_prefix\.sam 2>>$output_prefix\_bwa_stderr.log",$log);
		}
	} else {
		if ($options->{disable_SW}) {
			run_cmd("bwa sampe -s $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 2>>$output_prefix\_bwa_stderr.log | samtools view -bhS - > $output_prefix\.bam ",$log);
		} else {
			run_cmd("bwa sampe $ref $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 2>>$output_prefix\_bwa_stderr.log | samtools view -bhS - > $output_prefix\.bam ",$log);  ## 2>>$output_prefix\_bwa_stderr.log
		}
	}

  	## Sort and index bams 
  	if ($options->{sort_index}){
  		run_cmd("samtools sort $output_prefix\.bam $output_prefix\.srt",$log);
  		run_cmd("samtools index $output_prefix\.srt.bam $output_prefix\.srt.bai",$log);
  	}

    ## Samtools Mpileup
    if ($options->{mpileup}){
    	run_cmd("samtools mpileup -Af $options->{ref} $output_prefix\.srt.bam > $output_prefix\.COVERAGE.txt",$log);
    }

	## Picard insert metrics
	if ($options->{insert_metrics}){
		run_cmd("java -Xmx3g -jar /home/jdhotopp/bin/Picard/picard-tools-1.48/CollectInsertSizeMetrics.jar AS=true I=$output_prefix.srt.bam O=$output_prefix\_std_insert.metrics H=$output_prefix\_std_insert.histogram M=0 VALIDATION_STRINGENCY=SILENT",$log);
		run_cmd("java -Xmx3g -jar /home/jdhotopp/bin/Picard/picard-tools-1.48/CollectInsertSizeMetrics.jar AS=true I=$output_prefix.srt.bam O=$output_prefix\_lrg_insert.metrics H=$output_prefix\_lrg_insert.histogram M=0 VALIDATION_STRINGENCY=SILENT DEVIATIONS=1000000000000000000",$log);
	}

    ## Cleanup intermediate files (.sai)
    unless ($options->{no_cleanup}){
    	run_cmd("rm $output_prefix\.1.sai",$log);
    	run_cmd("rm $output_prefix\.2.sai",$log);
    	run_cmd("rm $output_prefix\_bwa_stderr.log",$log);
    }
    my $retval = defined $options->{sort_index} ? "$output_prefix\.srt.bam" : "$output_prefix\.bam";

    print STDERR "====== Completed BWA mapping: $file1 against: $ref output: $output_prefix\.bam ======\n";
    return $retval;
}

1