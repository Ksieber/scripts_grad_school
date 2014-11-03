#!/usr/local/bin/perl -w
use linecount;
use errorcheck;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;
my %options;
my $results = GetOptions (\%options,
         'bam=s',
         'bam_list=s',
         'not_sorted',
         'region=s',
         'region_file=s',
         'output_dir=s',
         'M_M',
         'UM_M',
         'discord',
         'wrong_dir',
         'sort_name',
         'help|h'
      );

if($options{help}){
   die "HELP: This script will parse a bam for the region(s) specified in the region_file.
      --bam= ~Required. 
      --bam_list=  ~Required. List of bams to pull the region(s) from. 
      --not_sorted Use this if the input bam(s) are NOT sorted by position already.
      --region=\tchr{x}:1-2
      --region_file= Name of the file with the desired regions to parse out.
      --M_M    Toggle on to find M_M reads in the region of interest.
      --UM_M   Toggle on to find reads UM_M to the region(s) of interest.
      --discord   Use this to pull out discordant reads within the desired regions.
      --wrong_dir Pull reads facing the wrong directions. (RR, FF, RF)
      --sort_name\tToggle on to sort the output bam by read name instead of by the default by position.
      --output_dir=  Specify where to put the output. Default = current directory.
     \n";
}

if(!$options{bam} && !$options{bam_list}){
   die "ERROR: Must give an input bam with --bam or --bam_list\n";
}
my $number_parsing_conditions = 0; 
if($options{M_M}){
   $number_parsing_conditions++;
}
if($options{UM_M}){
   $number_parsing_conditions++;
}
if($options{discord}){
   $number_parsing_conditions++;
}
if($options{wrong_dir}){
   $number_parsing_conditions++;
}
if($number_parsing_conditions>1){
   die "ERROR: Can't parse on multiple parameters at once. Use only one of the following options at once: M_M, UM_M, discordant, wrong_dir.\n";
}

my @desired_regions;
my @tmp_bam_list;
my @bam_list;
#################################################
## Setting up input and outputs
if($options{bam_list}){
   open(LIST, "<", "$options{bam_list}") or die "Unable to open input bam list: $options{bam_list} because: $!\n";
   while(<LIST>){ 
      chomp;
      push(@tmp_bam_list,$_);
   }
   close LIST or die;
}     
if($options{bam}){
   push(@tmp_bam_list,$options{bam});
}
##Sort the bams files if they are not sorted
if($options{not_sorted}){
    foreach my $bam (@tmp_bam_list){
        # $bam =~ m/(\w+?)\.bam$/;
        my($file,$dir,$suf)=fileparse($bam,".bam");
       `samtools sort $bam $dir/$file\_position_sorted`;
       &errchk($?);
       `samtools index $dir/$file\_position_sorted.bam $dir/$file\_position_sorted.bai`;
       &errchk($?);
       push(@bam_list,"$dir/$file\_position_sorted.bam");
    }
} else {
    @bam_list = @tmp_bam_list;
}
@tmp_bam_list = ();

if($options{region}){
   push(@desired_regions,$options{region});
}
if($options{region_file}){
   open(REGION, "<", "$options{region_file}") or die "Unable to open file with desired regions: $options{region_file} because: $!\n";
   while(<REGION>){
      chomp;
      push(@desired_regions, $_);
   }
   close REGION or die "Unable to close file: $options{region_file} because: $!\n";
}

#################################################
foreach my $bam (@bam_list){
   print STDERR "Processing: $bam . . .\n";
   my($output,$path,$suf)=fileparse($bam,".bam");
   my $outputdir = $options{output_dir} ? $options{output_dir} : $path;
   $outputdir=~s/\/^//;
   open(OUT, ">", "$outputdir/$output\_regions.reads");
   print STDERR "Making GOOD list with read names from regions of interest . . .\n";
   foreach my $region (@desired_regions){
      print STDERR "Parsing the input bam: $bam \tfor region: $region . . .\n";
      my $check=`samtools view $bam \'$region\' | wc -l`;
      chomp($check);
      print STDERR "Found: $check reads in the region: $region from the bam: $bam\n";
      if($options{UM_M}){
         `samtools view -f0x8 $bam \'$region\' | cut -f1 >> $outputdir/$output\_regions.reads`;
      } elsif ($options{M_M}){
        `samtools view -F0xC $bam \'$region\' | cut -f1 >> $outputdir/$output\_regions.reads`; 
      } elsif ($options{wrong_dir}){
         $old_wc=linecount::wc("$outputdir/$output\_regions.reads");
         `samtools view -f0x30 -F0xC $bam \'$region\' | cut -f1 >> $outputdir/$output\_regions.reads`;   ## RR
         my $new_wc=linecount::wc("$outputdir/$output\_regions.reads"); 
         $diff_wc=$new_wc-$old_wc;
         print STDERR "RR reads found: $diff_wc\n";
         $old_wc=$new_wc;
         `samtools view -uf0x30 $bam \'$region\' | samtools view -F0xC - | cut -f1 >> $outputdir/$output\_regions.reads`; ## FF
         $new_wc=linecount::wc("$outputdir/$output\_regions.reads");
         $diff_wc=$new_wc-$old_wc;
         print STDERR "FF reads found: $diff_wc\n";
         $old_wc=$new_wc;
         open(RF,"-|","samtools view -f0x10 -F0xC $bam \'$region\'");  ## Parse for RF
         while(<RF>){
            my @f=split; 
            if($f[8]>0){            ## if R read is upstream of F
               print OUT "$f[0]\n";
            }
         }
         close RF;
         $new_wc=linecount::wc("$outputdir/$output\_regions.reads");
         $diff_wc=$new_wc-$old_wc;
         print STDERR "RF reads found: $diff_wc\n";
         `sort $outputdir/$output\_regions.reads | uniq > $outputdir/tmp_regions.reads`;
         `mv $outputdir/tmp_regions.reads $outputdir/$output\_regions.reads`;
      }else {
         open(my $in, "-|", "samtools view $bam \'$region\'") or die "Unable to open file: $bam because: $!\n";
         while(<$in>){
            chomp;
            my @f=split;
            if($options{discord}){
               if($f[2] ne $f[6] && $f[6] !~/\=/){
                  print OUT "$f[0]\n";
               }
               if(abs($f[8]) > 10000){
                  print OUT "$f[0]\n";
               }
            }
            if(!$options{discord}){
               print OUT "$f[0]\n";
            }
         }
      }
   }
   print STDERR "\tdone making GOOD list.\nParsing the bam for the good reads\n";   
   `perl /home/ksieber/scripts/parse_bam_reads_from_list.pl --input=$bam --good_list=$outputdir/$output\_regions.reads --output_prefix=$output\_regions --output_dir=$outputdir`;
   if($options{sort_name}){
      `samtools sort -n $outputdir/$output\_regions.bam $outputdir/$output\_regions_sort_names`;
      &errchk($?);
      `rm $outputdir/$output\_regions.bam`;
    }
}
close OUT;
__END__
   















