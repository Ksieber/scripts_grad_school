#!/usr/local/bin/perl
use warnings;
use strict;
use File::Basename;
use setup_input;
use setup_output;
use errorcheck;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,                       
						'input=s',
                    	'input_list=s',
						'output_prefix=s',
                        'output_dir=s',
                        'subdirs',
                        'ref=s',
                        'ref_list=s',
                        't=s',
                        's',
                        'bam_output',
                        'sort_index_bams',
                        'mpileup',
                        'cleanup',
                        'insert_metrics',
                        'mapped_only',
                        'help|h',
);


if ($options{help}){die "\nHELP: This script will align the input (fastq/bam) to a reference.
\t--input=		Input file to be BWA mapped. Either: fastq, or bam. 
\t--input_list=		List of input files to be mapped. 1 bam/line. _1,_2 fastq/line (fastqs MUST be comma seperated).
\t--ref=			Reference.fna+index
\t--ref_list=		List of References.
\t--output_prefix=	Prefix for each output.  Ie. SRA_LGT_NC001234
\t--output_dir=		Will output to current working directory unless another is specified with this. ie. /ksieber_dir/tmp/
\t--subdirs		Turn this on to make subdirectories for each input file to be mapped.
\t--t=			# [1] Will set the number of cpu threads to use for bwa aln steps. USE CAREFULLY.   
\t--s			Turn this on to disable Smith-Waterman for the unmapped mate. Ideal for quicker LGT mappings IF they are high confidence. 
\t--mapped_only		Keep only mates with 1 mapped read.
\t--bam_output		Convert the new.sam into a new.bam
\t--sort_index_bams	Sort and index the new.bam into new.srt.bam and new.srt.bai.  If --bam_output wasn't used it will generate the .bam anways. 
\t--mpileup		Calculate pileup coverage on .bam.  Must be used with --bam_output
\t--cleanup		WARNING: This option delete new.1.sai; new.2.sai; Leaving only new.srt.bam/bai
\t--insert_metrics	Use Picard to calculate insert size metrics. Must be used with --bam_output.
\t--help\n";
}

if(!$options{input} || !$options{input_list}){die "Error: Must give input files to map with --input or --input_list.\n";}
if (!$options{ref} && !$options{ref_list}){die "ERROR:  Must have enter a reference file to use.\n";}
if ($options{insert_metrics} && !$options{bam_output}){die "ERROR: Must use --bam_output with --insert_metrics.\n";}
my @in_suffix_list=('.bam','.fastq','.fq');  
my @ref_suffix_list=('.fa','.fna','.txt');
my $threads = $options{t} ? $options{t} : "1";
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
my $out_dir=setup_output($input);

## Run BWA
foreach my $files (@$input){
	foreach my $refs (@ref_list){
		bwa_align{$files,$refs};
	}
}

sub bwa_align {
  	my ($files,$ref) = @_;
	my ($ref_name,$ref_path,$ref_suf)=fileparse($ref,@ref_suffix_list); 
	my $file1;
	my $file2;
	my $bam;     															# 0= fastq, 1=bam
	if($files=~/.fq$/ || $files=~/.fastq$/){
		($file1,$file2)=split(/,/,$files);
		$bam=0;
	} elsif ($files=~/.bam$/){
		$file1=$files;
		$file2=$files;
		$bam=1;
	} else {
		die "Could not resolve the input type. Make sure it is either a .bam,.fq,.fastq\n";
	}
	my ($input,$path,$suf)=fileparse($file1,@in_suffix_list);
	my $out = $options{output_prefix} ? "$options{output_prefix}\_$ref" : "$input\_at_$ref";
	my $dir = $out_dir->{$files}; 
  	$dir =~s/\/$//g;
  	my $output_prefix= "$dir\/$out";
  	if ($bam==1){
    	`bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref -b1 $file1 > $output_prefix\.1.sai`;
	   	 &errchk($?);
    	`bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref -b2 $file2 > $output_prefix\.2.sai`;
    	&errchk($?);
  	} else {
    	`bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref $file1 > $output_prefix\.1.sai`;
	    &errchk($?);
    	`bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $ref $file2 > $output_prefix\.2.sai`;
	    &errchk($?);
	  }
	if ($options{mapped_only}){
    	`bwa sampe $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 | samtools view -F0x4 -bhS - > $output_prefix\.bam`;
  	} elsif ($options{bam_output} or $options{sort_index_bams}){
	    if($options{s}){
    		`bwa sampe -s $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 | samtools view -bhS - > $output_prefix\.bam`;
    	} else {
    		`bwa sampe $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 | samtools view -bhS - > $output_prefix\.bam`;
    	}
  	} else {
    	if($options{s}){
    		`bwa sampe -s $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 > $output_prefix\.sam`;
    	} else {
    		`bwa sampe $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $file1 $file2 > $output_prefix\.sam`;
    	}
  	}
  	&errchk($?);
  	if ($options{sort_index_bams}){
    	`samtools sort $output_prefix\.bam $output_prefix\.srt`;
    	&errchk($?);
    	`samtools index $output_prefix\.srt.bam $output_prefix\.srt.bai`;
    	&errchk($?);
    }
    if ($options{bam_output} && $options{mpileup}){
    	`samtools mpileup -Af $options{ref} $output_prefix\.srt.bam > $output_prefix\.COVERAGE.txt`;
      	&errchk($?);
	}
    if ($options{insert_metrics}){
 		`java -Xmx3g -jar /home/jdhotopp/bin/Picard/picard-tools-1.48/CollectInsertSizeMetrics.jar AS=true I=$output_prefix.srt.bam O=$output_prefix\_std_insert.metrics H=$output_prefix\_std_insert.histogram M=0 VALIDATION_STRINGENCY=SILENT`;
      	&errchk($?);
		`java -Xmx3g -jar /home/jdhotopp/bin/Picard/picard-tools-1.48/CollectInsertSizeMetrics.jar AS=true I=$output_prefix.srt.bam O=$output_prefix\_lrg_insert.metrics H=$output_prefix\_lrg_insert.histogram M=0 VALIDATION_STRINGENCY=SILENT DEVIATIONS=1000000000000000000`;
		&errchk($?);
    }
	if ($options{cleanup}){
    	`rm $output_prefix\.1.sai`;
	    `rm $output_prefix\.2.sai`;
	    if ($options{sort_index_bams}){
    	  `rm $output_prefix\.bam`;
    	}
  	}
}

__END__
