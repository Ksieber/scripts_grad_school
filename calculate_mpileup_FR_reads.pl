#!/usr/local/bin/perl
#use strict;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
      'input=s',
      'input_list=s',
      'region=s',
      'output_prefix=s',
      'output_dir=s',
      'M_M',
      'M_UM',
      'A=s',
      'discord',
      'wrong_direction',
      'help',
      );

if($options{help}){die
   "HELP: This script will take a bam file and calculate the coverage of the Foward and Reverse reads.
      --input=          Input bam
      --input_list=     List of bams to process.
      --region=         chr#:100-200. Use to look @ reads only in this region.
      --output_prefix=  [Optional] Prefix for the output file. Suffix = .F_mpileup && .R_mpileup
      --output_dir=     [Optional] Name the directory for output. Default = cwd
      --M_M             [Optional] Toggle on to only calculate Mapped_Mapped reads only. Highly recommended to use. 
      --M_UM            [Optional] Toggle on to only calculate Mapped_UN-Mapped reads only. Works for LGT reads.
      --wrong_direction [Optional] Toggle on to only calculate F-F, R-R, & R-F 
      --A=              =<0|1> Override samtools mpileup -A. Default = 1 (ON).\n";
}
#### Check input
if(!$options{input} && !$options{input_list}){die "ERROR: Must give an input bam with --input=<BAM>. Try again.\n";}
if($options{M_M} && $options{M_UM}){die "ERROR: Can only use one: M_M or M_UM at once. Try again.\n";}


### Make a list of input bams
my @in_list;
if($options{input_list}){
   open(POO,"<","$options{input_list}") || die "ERROR: Couldn't open the input list: $options{input_list} because: $!\n";
   while(<POO>){
      chomp;
      push(@in_list,$_);
   }
   close POO;
}
if($options{input}){
   push(@in_list,$options{input});
}

### Process Each bam
foreach my $bam (@in_list){
   ##### Assign output names
   my ($filename,$directory) = fileparse($bam, ".bam");
   my $prefix = $options{output_prefix} ? $options{output_prefix} : $filename;
   my $dir = $options{output_dir} ? $options{output_dir} : $directory;
   my $out = "$dir/$prefix";

   #### Determine options for samtools view & mpileup
   my $view = "-hu";              ## Default
   my $mpileup = "-Ad 100000";    ## Default
   my $f = "-F0x10";              ## samtools view options for looking @ Forward stranded reads
   my $r = "-f0x10";              ## samtools view options for looking @ Reverse stranded reads
   my $region = $options{region} ? $options{region} : " " ;
   
   if($options{wrong_direction}){      ## Had to make a seperate section for processing because we need to do serial parsing for FF &| RR. 
      $view = "-huF0xC";
      $mpileup = "-Ad 100000";
      my $ff = "-F0x30";
      my $rr = "-f0x30";
      my $check = 0;
      $check = `samtools view $ff $bam $region | wc -l`;
      chomp($check);
#print STDERR "samtools view $ff $bam $region | wc -l\n";
#print STDERR "CHECK1: $check\n"; 
      if($check!=0){
         `samtools view $view $ff $bam $region | samtools mpileup $mpileup - > $out.FF_mpileup`;     ## looking @ FF reads
      } else {print STDERR "No FF reads found for $bam\n";}
      $check = 0;
      $check = `samtools view $rr $bam $region | wc -l`;
      chomp($check);
#print STDERR "CHECK2=$check\n"; 
      if($check!=0){
         `samtools view $view $rr $bam $region | samtools mpileup $mpileup - > $out.RR_mpileup`;     ## looking @ RR reads
      } else {print STDERR "No RR reads found for $bam\n";}
      $check = 0;
      $check = `samtools view $view $r $bam $region | samtools view -F0x20 - | wc -l`;
      chomp($check);
#print STDERR "CHECK3=$check\n"; 
      if($check!=0){
         `samtools view $view $r $bam $region | samtools view -hF0x20 - | /home/ksieber/scripts/RF_mpileup.sh | samtools view -Shu - | samtools mpileup $mpileup - > $out.RF_mpileup`;
      } else {print STDERR "No RF reads foudn for $bam\n";}
      next; 
   }
   if($options{M_M}){
      $view = "-huF0xC";
      $mpileup = "-d 100000";
   }
   if($options{M_UM}){
      $view = "-huf0x8 -F0x4";
      $mpileup = "-Ad 100000";
   }  
   if($options{A} eq 0){
      $mpileup = "-d 100000";
   }
   `samtools view $view $f $bam $region | samtools mpileup $mpileup - > $out.F_mpileup`;
   `samtools view $view $r $bam \'$region\' | samtools mpileup $mpileup - > $out.R_mpileup`;
}

sub wc { 
   $_[0]=~/^(\d+)\s+\w+/;
   my $num=$1;
   return "$num";
}
