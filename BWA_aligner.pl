#!/usr/local/bin/perl
use strict;
use File::Basename;
use lib '/home/ksieber/scripts/';
use run_cmd;
use setup_input;
use setup_output;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
our $results = GetOptions (\%options,                       
						'input=s',
                    	'input_list=s',
						'output_prefix=s',
                        'output_dir=s',
                        'subdirs=s',
                        'ref=s',
                        'ref_list=s',
                        't=s',
                        'disable_SW=s',
                        'bam_output=s',
                        'sort_index_bams=s',
                        'mpileup=s',
                        'no_cleanup=s',
                        'insert_metrics=s',
                        'mapped_only=s',
                        'cmd_log=s',
                        'help|h',
);

if ($options{help}) {die "\nHELP: This script will align the input (fastq/bam) to a reference.
\t--input=			Input file to be BWA mapped. Either: fastq, or bam. (Fastq_1,Fastq_2)
\t--input_list=			List of input files to be mapped. 1 bam/line. _1,_2 fastq/line (fastqs MUST be comma seperated).
\t--ref=				Reference.fna+index
\t--ref_list=			List of References.
\t--output_prefix=		Prefix for each output.  Ie. (SRA_LGT)_at_\$ref_name
\t--output_dir=			Will output to current working directory unless another is specified with this. ie. /ksieber_dir/tmp/
\t--subdirs=			<0|1> [0] 1=Make subdirectories for each input file to be mapped.
\t--t=				<#>   [1] Set the number of cpu threads to use for bwa aln steps. USE CAREFULLY.   
\t--disable_SW=			<0|1> [0] 1=Disable Smith-Waterman for the UM mate. Ideal for quicker LGT mappings IF they are high confidence. 
\t--mapped_only=			<0|1> [0] 1=Only keep mates with 1 mapped read.
\t--sam_output=			<0|1> [0] 0=.bam output; 1=.sam output
\t--sort_index_bams=		<0|1> [0] 1=Sort and index the new.bam into new.srt.bam and new.srt.bai. 
\t--mpileup=			<0|1> [0] 1=Calculate pileup coverage on .bam.
\t--no_cleanup=			<0|1> [0] 0=Removes .sai files and unsorted.bam with --sort_index_bams. 1=No deleting intermediate data. 
\t--insert_metrics=		<0|1> [0] 1= Use Picard to calculate insert size metrics.
\t--cmd_log=			<0|1> [0] 1= Log all commands run in each output_dir/output_prefix.cmd_log
\t--help\n";
}

if (!$options{input} && !$options{input_list}) {die "Error: Must give input files to map with --input or --input_list.\n";}
if (!$options{ref} && !$options{ref_list}) {die "ERROR:  Must have enter a reference file to use.\n";}
my @in_suffix_list=('.bam','.fastq','.fq');  
my @ref_suffix_list=('.fa','.fna','.txt');
my $threads = $options{t} ? $options{t} : "1";													## Default # of threads = 1
if ($options{mpileup}==1) {$options{sort_index_bams}=1;}										## Mpileup needs a srt.bam, so turn it on automatically if --mpileup is passed
my @ref_list;


## Setup the reference list
if($options{ref}){
	push(@ref_list,$options{ref});
}
if($options{ref_list}){
	open(LIST,"<","$options{ref_list}") || die "Error: Can't open reference list because: $!\n";
	while(<LIST>){
		chomp;
		push(@ref_list,$_);
	}
	close LIST;
}

## Setup the input list
my $input=setup_input();

## Setup output directories
my $out_dir;
for my $first_input (@$input[0]){
	if($first_input=~/\.bam$/){
		$out_dir=setup_output($input);
	}
	if($first_input=~/\.fq$/ || $first_input=~/\.fastq$/){
		my @fastq_list;
		foreach my $fastqs (@$input){
			my ($file1,$file2)=split(/,/,$fastqs);
			push(@fastq_list,$file1);
		}
		$out_dir=setup_output(\@fastq_list);
	}
}	

## Setup Command logs
my $cmd_logs; 
if ($options{cmd_log}==1) {
	$cmd_logs=setup_logs($out_dir);
}
  		
## Run BWA
foreach my $files (@$input){
	foreach my $refs (@ref_list){
		bwa_align($files,$refs);
	}
}

## BWA alignment
sub bwa_align {
  	my ($files,$ref) = @_;
	my ($ref_name,$ref_path,$ref_suf)=fileparse($ref,@ref_suffix_list);     ## Grab the reference name to use for the naming the output
	my $file1;																## Global input file name
	my $file2;																## Global input file name2
	my $bam;     															## 0=fastq, 1=bam
	if($files=~/.fq$/ || $files=~/.fastq$/){
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
	if($options{cmd_log}==1){
		$log=$cmd_logs->{$file1};
	}
	
	## Setup output prefix (path/file-name).bam
	my ($input,$path,$suf)=fileparse($file1,@in_suffix_list);
	my $out = $options{output_prefix} ? "$options{output_prefix}\_$ref_name" : "$input\_at_$ref_name";
	my $dir = $out_dir->{$file1}; 
	if ($dir=~/\/$/) {$dir =~s/\/$//g;}
  	my $output_prefix= "$dir\/$out";
  	
  	
  	## Run BWA ALN
  	if ($bam==1) {
    	run_cmd("bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref -b1 $file1 > $output_prefix\.1.sai",$log);
    	run_cmd("bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref -b2 $file2 > $output_prefix\.2.sai",$log);
  	} else {
    	run_cmd("bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref $file1 > $output_prefix\.1.sai",$log);
    	run_cmd("bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref $file2 > $output_prefix\.2.sai",$log);
	}
	
	## BWA SAMPE
	if ($options{mapped_only}==1){
    	run_cmd("bwa sampe $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 | samtools view -F0x4 -bhS - > $output_prefix\.bam",$log);
  	} elsif ($options{sam_output}==1) {
    	if ($options{disable_SW}==1) {
    		run_cmd("bwa sampe -s $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 > $output_prefix\.sam",$log);
    	} else {
    		run_cmd("bwa sampe $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 > $output_prefix\.sam",$log);
    	}
  	} else {
	    if ($options{disable_SW}==1) {
    		run_cmd("bwa sampe -s $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 | samtools view -bhS - > $output_prefix\.bam",$log);
    	} else {
    		run_cmd("bwa sampe $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 | samtools view -bhS - > $output_prefix\.bam",$log);
    	}
  	}
  	
  	## Sort and index bams 
  	if ($options{sort_index_bams}==1){
    	run_cmd("samtools sort $output_prefix\.bam $output_prefix\.srt",$log);
    	run_cmd("samtools index $output_prefix\.srt.bam $output_prefix\.srt.bai",$log);
    }
    
    ## Samtools Mpileup
    if ($options{mpileup}==1){
    	run_cmd("samtools mpileup -Af $options{ref} $output_prefix\.srt.bam > $output_prefix\.COVERAGE.txt",$log);
	}
	
	## Picard insert metrics
    if ($options{insert_metrics}==1){
 		run_cmd("java -Xmx3g -jar /home/jdhotopp/bin/Picard/picard-tools-1.48/CollectInsertSizeMetrics.jar AS=true I=$output_prefix.srt.bam O=$output_prefix\_std_insert.metrics H=$output_prefix\_std_insert.histogram M=0 VALIDATION_STRINGENCY=SILENT",$log);
		run_cmd("java -Xmx3g -jar /home/jdhotopp/bin/Picard/picard-tools-1.48/CollectInsertSizeMetrics.jar AS=true I=$output_prefix.srt.bam O=$output_prefix\_lrg_insert.metrics H=$output_prefix\_lrg_insert.histogram M=0 VALIDATION_STRINGENCY=SILENT DEVIATIONS=1000000000000000000",$log);
    }
    
    ## Cleanup intermediate files (.sai)
	unless ($options{no_cleanup}==1){
    	run_cmd("rm $output_prefix\.1.sai",$log);
	    run_cmd("rm $output_prefix\.2.sai",$log);
	    if ($options{sort_index_bams}==1){
    	  run_cmd("rm $output_prefix\.bam",$log);
    	}
  	}
}

__END__
