#!/usr/local/bin/perl5 -w
use strict;
use File::Basename;
if(@ARGV!=1){die "Must pass ARGV an input fasta file to substitute the id's for.\n";}
my ($file,$dir,$suf) = fileparse($ARGV[0],".fa");
$dir=~s/\/$//;
if($dir!~/\w+/){die "Must give full file path to input.\n";}
`perl /home/ksieber/scripts/sub_in_fa_unique_ids.pl $ARGV[0]`;
`clustalw -INFILE=$dir/$file\_sub-id.fa -OUTPUT=PHYLIP -OUTFILE=$dir/$file\_MSA.phy`;
my $random = int(rand(1000)) + 1;
open(OUT,">","Phyml_$random\.sh");
print OUT "phyml -i $dir/$file\_MSA.phy -b 1000\n";
close OUT;
`qsub -P jdhotopp-lab Phyml_$random\.sh`;

