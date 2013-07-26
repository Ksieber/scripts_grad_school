#!/usr/local/bin/perl -w
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

## This script will take SEPERATE Bacterial & Hg19 LGT.sam files can convert them into a Links.txt for circos
## This script also needs a tab-list of bacterial chromosomes and the links name for each chr, Ex:
## Bac1   NC_1       (no .1)
## Bac2   NC_2

my %options;
my $results = GetOptions (\%options,
			  'bac_sam=s',
			  'hg_sam=s',
			  'host_sam=s',
			  'chr_list=s',
			  'output=s',    ##specify the output name instead of parsing it out automatically, useful for designating a different directory for output      
			  'help|h',
			 );
##  'alt_output',  ##use for circos read coverage histogram; reverse hg-->bac links instead of normal output (bac-->hg)

if ($options{help}){die "\nHELP: This script will parse two LGT.sam files (Hg19 & Bacterial) and create the links file needed for Circos.
   This script also needs a tab delimited list of the chr's from the bacteria and the names to call them in the links file.  EX:
      Bac1  NC_1
      Bac2  NC_2
      Lambda\tNC_001416
   
   --bac_sam *            (Mandatory)
   --hg_sam  *            (Mandatory)
   --host_sam            (~Mandatory)  Host for transfers.  ie if not hg19, could be Dana, HeLa etc
   --chr_list*            (Mandatory)
   --output               (optional, override output name)
   --alt_output           (optional, used for circos read coverage histogram)
   --help\n";}


## Varialbes
my $output;
my %bac_counter;
my %bac_output;
my $bac_chr;
my %hg_counter;
my %hg_output;
my $read_length;
my $second_position;
my $hgChr;
my %chr_converter;
my $chr_sam_name;
my $chr_links_name;
my $alt_output;


## Name output based on input file
if ($options{bac_sam}){
    $options{bac_sam} =~ m/^(.*)\.sam/;
    $output = $options{output} ? $options{output} : "$1\.Links.txt";
    if ($options{alt_output}){
	$alt_output = "$output.ALT";
    }
}
 
if ($options{mixed_sam}){
    $options{mixed_sam} =~ m/^(.*)\.sam/;
    $output = $options{output} ? $options{output} : "$1\.Links.txt";
    if ($options{alt_output}){
	$alt_output = "$output.ALT";
    }
}
 

## Make hash to convert chr_sam names into links_sam_names
if (!$options{chr_list}){die "ERROR: Must give chr_list\n";}
open (LIST, "<", $options{chr_list}) || die "ERROR: Couldn't open $options{input_sam} because: $!\n";
while (<LIST>){
  ($chr_links_name,$chr_sam_name)=(split)[0,1];
  $chr_converter{$chr_sam_name}=$chr_links_name;
}
close LIST || die "ERROR: Couldn't close $options{chr_list} because: $!\n";

open (OUTPUT, ">", $output);
if ($options{alt_output}){
   open (ALTOUTPUT, ">", $alt_output);
}

if ($options{bac_sam}){
  open (BAC_SAM, "<", $options{bac_sam}) || die "Can't open $options{bac_sam} because: $!\n";
  while (<BAC_SAM>){
    my ($header,$chr,$position,$bases)=(split)[0,2,3,9];
    next if ($header=~m/^\@SQ/);
    if ($chr =~ m/gi\|.*\|(NC\_.*)\..\|*/ && ($chr_converter{$1})){               ## If a bacaterial read
      $bac_counter{$header}++;
      if ($bac_counter{$header}==1){               ## process first bacterial read
	$read_length = length($bases);             ## Calculate read length
	$second_position = $position + $read_length;  
	$bac_output{$header} = "$header\t$chr_converter{$1}\t$position\t$second_position\n";     ## Save line for later
      } 
    }
  }
  close BAC_SAM || die "ERROR: Couldn't close $options{bac_sam} because $!\n";
}

if ($options{hg_sam}){
  open (HG_SAM, "<", $options{hg_sam}) || die "ERROR: Couldn't open $options{hg_sam} because: $!\n";
  while(<HG_SAM>){
    my ($header,$chr,$position,$bases)=(split)[0,2,3,9];
    next if ($header=~m/^\@SQ/);
    if ($chr =~ m/chr(.{1,2})$/){                  ## If Human 1-22+x+y    ## Remember the chr name
      $hg_counter{$header}++;
      if ($hg_counter{$header}==1){                ## Process the first read
	      $read_length = length($bases);
      	$second_position = $position + $read_length;
	      $hg_output{$header} = "$header\ths$1\t$position\t$second_position\n";
      }
    } 
    if ($bac_output{$header} &&  $hg_output{$header} && $bac_counter{$header} >=2 && $hg_counter{$header} >= 2) {
      print OUTPUT "$bac_output{$header}$hg_output{$header}";
      if ($options{alt_output}){
	print ALTOUTPUT "$hg_output{$header}$bac_output{$header}";
      }
    } 
  }
  close HG_SAM || die "ERROR closing $options{input_sam}\n";
  close OUTPUT || die "ERROR closing $output\n";
  if ($options{alt_output}){
    close ALTOUTPUT || die "ERROR closing $alt_output\n";
  } 
}

if ($options{host_sam}){
  open (HOST_SAM, "<", $options{host_sam}) || die "ERROR: Couldn't open $options{host_sam} because: $!\n";
  while(<HOST_SAM>){
    my ($header,$chr,$position,$bases)=(split)[0,2,3,9];
    if ($chr =~ m/gi\|.*\|..*\|(.*)\..\|$/){                 ## Remember the chr name
      next if ($header=~m/^\@SQ/);
      next if (!$chr_converter{$1});
      $hg_counter{$header}++;
      if ($hg_counter{$header}==1){                ## Process the first read
	$read_length = length($bases);
	$second_position = $position + $read_length;
	$hg_output{$header} = "$header\t$chr_converter{$1}\t$position\t$second_position\n";
      }
    } 
    if ($bac_output{$header} &&  $hg_output{$header} && $bac_counter{$header} >=2 && $hg_counter{$header} >= 2) {
      print OUTPUT "$bac_output{$header}$hg_output{$header}";
      if ($options{alt_output}){
	print ALTOUTPUT "$hg_output{$header}$bac_output{$header}";
      }
    } 
  }
  close HOST_SAM || die "ERROR closing $options{host_sam}\n";
  close OUTPUT || die "ERROR closing $output\n";
}
#if ($options{mixed_sam}){
 #   open (CLONE, "<", $options{mixed_sam}) || die "ERROR: Couldn't open $options{by_clone} for reading b/c: $!\n";
  #  while(<CLONE>){
#	my ($header,$chr,$position,$cigar,$bases)=(split)[0,2,3,5,9];
#	next if ($header=~m/^\@SQ/); ## && print STDERR "Fail: Began w/ \@SQ\n";
#	next if ($cigar !~ m/M/); ## && print STDERR "Fail: Didn't Map\n";
#	$chr =~ m/gi\|.*\|..*\|(.*)\..\|$/;
#	next if (!$chr_converter{$1}); ## && print STDERR "Fail: $1 not on the chr.list\n";
#	##print STDERR "test2\t$1\t$chr_converter{$1}\n";
#	$read_length = length($bases);             ## Calculate read length
#	$second_position = $position + $read_length;  
#	print OUTPUT "$header\t$chr_converter{$1}\t$position\t$second_position\n";    
 #   } 
#}
#close CLONE;




__END__
    
