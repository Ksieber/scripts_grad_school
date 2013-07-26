#!/usr/local/bin/perl
use warnings;
use strict;
##This script will doing bwa mapping
##This script will also use samtools to create the .bam files
use File::Basename;
use strict;
use errorcheck;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
                          'fastq_1=s',
                          'fastq_2=s',
                          'output_prefix=s',
                          'bam_input=s',
                          'bam_list=s',
                          'fastq_list=s',
                          'output_dir=s',
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
\t--output_prefix=     Prefix for each output.  Ie. SRA_LGT_NC001234
\t--output_dir=        Will output to current working directory unless another is specified with this. ie. /ksieber_dir/tmp/
\t--ref=               Reference.fna+index
\t--ref_list=          List of References.
\t--fastq_1=	        fastq_1 file
\t--fastq_2=	        fastq_2 file
\t--bam_input=         data.bam
\t--bam_list=          List of input bams. 1 input.bam per line. 
\t--fastq_list=        List of fastq files to align.  1 pair of _1&_2 per line. Each line is a different fastq pair.        
\t--t=                 # [1] Will set the number of cpu threads to use for bwa aln steps. USE CAREFULLY.   
\t--s                  Disable Smith-Waterman for the unmapped mate. Ideal for quicker LGT mappings IF they are high confidence. 
\t--mapped_only        Keep only mates with 1 mapped read.
\t--bam_output         Convert the new.sam into a new.bam
\t--sort_index_bams    Sort and index the new.bam into new.srt.bam and new.srt.bai.  If --bam_output wasn't used it will generate the .bam anways. 
\t--mpileup            Calculate pileup coverage on .bam.  Must be used with --bam_output
\t--cleanup            WARNING: This option delete new.1.sai; new.2.sai; Leaving only new.srt.bam/bai
\t--insert_metrics     Use Picard to calculate insert size metrics. Must be used with --bam_output.
\t--help
";}

if (!$options{fastq_1} && !$options{fastq_2} && !$options{fastq_list} && !$options{bam_input} && !$options{bam_list}){die "ERROR: Must give an input to align.  Use either --fastq_1/_2, fastq_list, bam_input, or bam_list.\n";}
if (!$options{ref} && !$options{ref_list}){die "ERROR:  Must have enter a reference file to use.\n";}
if ($options{insert_metrics} && !$options{bam_output}){die "ERROR: Must use --bam_output with --insert_metrics.\n";}
my @in_suffix_list=('.bam','.fastq','.fq');  
my @ref_suffix_list=('.fa','.fna','.txt');
my $threads = $options{t} ? "$options{t}" : "1";
my $NC=0;
my $ref;

if ($options{ref_list}){                                    ## Map at multiple refs
  open (LIST, "<", $options{ref_list}) || die "ERROR: Can't open $options{ref_list} because: $!\n";
  my @ref_list;
  while(<LIST>){
    chomp;
    push(@ref_list, $_);
  }
  close LIST || die "ERROR: Can't close $options{ref_list} because: $!\n";
  foreach my $reference (@ref_list){
    my ($ref_name,$ref_path,$ref_suf)=fileparse($reference,@ref_suffix_list);      
    $ref = $ref_name;
    $options{ref} = $reference;
    if ($options{fastq_1} && $options{fastq_2}){
      &bwa_align($options{fastq_1},$options{fastq_2});
    }
    if ($options{bam_input}){
      &bwa_align($options{bam_input},$options{bam_input});
    }
    if($options{bam_list}){
      my @bam_list;
      open (INLIST, "<", $options{bam_list}) || die "ERROR: Can't open $options{bam_list} because: $!\n";
      while(<INLIST>){
	      chomp;
	      push(@bam_list, $_);
      }
      close INLIST || die "ERROR: Can't close $options{bam_list} because: $!\n"; 
      foreach my $bam_input (@bam_list){
	      $options{bam_input}=$bam_input;
	      &bwa_align($options{bam_input},$options{bam_input});
      }
    }
  }
} elsif ($options{fastq_1} && $options{fastq_2}){             ## Map 2 Fastq Files @ A reference
  &bwa_align($options{fastq_1},$options{fastq_2});
} elsif ($options{bam_input}){                                ## Map A bam input file @ A reference
   my ($ref_name,$ref_path,$ref_suf)=fileparse($options{ref},@ref_suffix_list);
   $ref = $ref_name;
   &bwa_align($options{bam_input},$options{bam_input});
} elsif ($options{bam_list}){                                 ## Map a list of bams @ A reference
  my @bam_list;
  open (INLIST, "<", $options{bam_list}) || die "ERROR: Can't open $options{bam_list} because: $!\n";
  while(<INLIST>){
    chomp;
    push(@bam_list, $_);
  }
  close INLIST || die "ERROR: Can't close $options{bam_list} because: $!\n";
  my ($ref_name,$ref_path,$ref_suf)=fileparse($options{ref},@ref_suffix_list); 
  $ref = $ref_name;
  foreach my $bam_input (@bam_list){
    $options{bam_input}=$bam_input;
    &bwa_align($options{bam_input},$options{bam_input});
  }
}


sub bwa_align {
  my ($fastq_1,$fastq_2) = ($_[0],$_[1]);
  my ($input,$path,$suf)=fileparse($fastq_1,@in_suffix_list);
  my $out = $options{output_prefix} ? "$options{output_prefix}\_at_$ref" : "$input\_at_$ref";
  my $dir = $options{output_dir} ? $options{output_dir} : $path; 
  $dir =~s/\/$//g;
  my $output_prefix= "$dir\/$out";
  if ($options{bam_input}){
    `bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $options{ref} -b1 $fastq_1 > $output_prefix\.1.sai`;
    &errchk($?);
    `bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $options{ref} -b2 $fastq_2 > $output_prefix\.2.sai`;
    if ($? != 0 ){die "Fatal Error line 124.\n";}
  } else {
    `bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $options{ref} $fastq_1 > $output_prefix\.1.sai`;
    &errchk($?);
    `bwa aln -e -1 -M 3 -E 4 -O 11 -t $threads -o 1 $options{ref} $fastq_2 > $output_prefix\.2.sai`;
    &errchk($?);
  }
  if ($options{mapped_only}){
    `bwa sampe $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $fastq_1 $fastq_2 | samtools view -F0x4 -bhS - > $output_prefix\.bam`
  } elsif ($options{bam_output} or $options{sort_index_bams}){
    if($options{s}){`bwa sampe -s $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $fastq_1 $fastq_2 | samtools view -bhS - > $output_prefix\.bam`;}
    else{`bwa sampe $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $fastq_1 $fastq_2 | samtools view -bhS - > $output_prefix\.bam`;}
  } else {
    if($options{s}){`bwa sampe -s $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $fastq_1 $fastq_2 > $output_prefix\.sam`;}
    else{`bwa sampe $options{ref} $output_prefix\.1.sai $output_prefix\.2.sai $fastq_1 $fastq_2 > $output_prefix\.sam`;}
  }
  &errchk($?);
  if ($options{sort_index_bams}){
    `samtools sort $output_prefix\.bam $output_prefix\.srt`;
    &errchk($?);
    `samtools index $output_prefix\.srt.bam $output_prefix\.srt.bai`;
    &errchk($?);
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
