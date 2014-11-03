#!/usr/local/bin/perl
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use IntervalTree;

my %options;
my $results = GetOptions (\%options,
                          'dfam_dir=s',
                          'dfam_file=s',
                          'tree_list=s',
                          'tree_dir=s',
                          'bam=s',
                          'bam_list=s',
                          'not_position_sorted',
                          'output_dir=s',
                          'help',
);

if($options{help}){die 
   "HELP: This script will take a bam file and search the Dfam hmm Transposable Element hits \"database\" to see if the read maps on a TE. 
\t--bam=\t\t\t<FILE> ~REQUIRED. This bam is sorted by position, mapped at hg19 and will be searched for TE hits. Assuming these are only M_M reads.   MIGHT NOT HAVE TO BE SORTED IF I CAN READ ALL DFAM INTO MEM
\t--bam_list=\t\t<LIST> ~Required. List of bam inputs, 1 per line. 
\t--not_position_sorted\tToggle this on if the bams are not position sorted. This will sort the bams by position. 
\t--output_dir=\t\tOptional. Directory to put the output. Defaults to current working directry.
\t--dfam_dir=\t\tOptional. Directory containing the Dfam hmm hits split by chromosome.
\t--dfam_file=\t\tOptional. This file contains the Dfam hmm hits as one entire file. NOT split by chr.
\t--tree_list=\t\tOptional. List of trees for Interval_Tree saved to disk. This will skip the build step and just read the tree into mem.
\t--tree_dir=\t\tOptional. Directory with all the tree files split by chr.\n";
}

my %current_chr;
my $dfam_dir = $options{dfam_dir} ? $options{dfam_dir} : "/local/projects/HLGT/ksieber_dir/Dfam_hits";
my $dfam_file = $options{dfam_file} ? $options{dfam_file} : "/local/projects/HLGT/ksieber_dir/Dfam_hits/Dfam.hits";
my $tree;
my @tmp_bam_list;
my @bam_list;

if(!$options{bam} && !$options{bam_list}){
    die "Must give an input to search for TE hits. Use --bam=<FILE> or --bam_list=<LIST> \n";
}

if($options{bam}){
    push(@tmp_bam_list,$options{bam});
}

if($options{bam_list}){
    open(LIST, "<", "$options{bam_list}") or die "Unable to open input bam list: $options{bam_list} becase $!\n";
    while(<LIST>){
        chomp;
        push(@tmp_bam_list,$_);
    }
    close LIST or die;
}

##Sort the bams files by position
if($options{not_position_sorted}){
    foreach my $bam (@tmp_bam_list){
        $bam =~ m/(\w+?)\.bam$/;
        `samtools sort $bam $1_position_sorted`;
        push(@bam_list,"$1_position_sorted.bam");
    }
} else {
    @bam_list = @tmp_bam_list;
}
@tmp_bam_list = ();


##########################################################################

## Process with prebuilt dfam trees split by chr.
if($options{tree_dir}){
    foreach my $bam (@bam_list){
        open(my $in, "-|", "samtools view $bam") or die "Unable to open the input bam file: $bam because $!\n";    ## Read in input BAM
        $bam =~ m/(\w+?)\.bam$/;
        open(OUT, ">", "$options{output_dir}$1.dfam_overlap") or die "Unable to open output $options{output_dir}/$1.dfam_overlap because $!\n";
        while(<$in>){
            chomp;
            my @line=split(/\t/,$_);                                                           ## bam line split
            #print STDERR "$line[2]\n";
            if($line[2] =~/chrM/){
               print OUT "$line[0]\tchrM\n";
               next;
            }
            if(!$current_chr{$line[2]}){                                                       ##If I haven't seen this chr yet
                print STDERR "loading tree: $options{tree_dir}/$line[2].tree . . .\n";
                $tree = IntervalTree->new();                                                   ##clear tree
                $tree->treeFromFile("$options{tree_dir}/$line[2].tree");                                                    ## 
                print STDERR "Tree loaded.\n";
                $current_chr{$line[2]}++;
            }
            if($current_chr{$line[2]}){                                                       ## If we already have the tree built for this chr
                print OUT "$line[0]\t";
                my @overlaps = $tree->searchInterval($line[3], $line[3]+length($line[9]));                          ## Find the intervals for each read
                if(@overlaps){                                                                ## If output, print it
                    foreach my $p (@overlaps){
                        print OUT "$p->[2]\t";
                    }
                } elsif (!$tree->searchInterval($line[3], $line[3]+length($line[9]))){
                    print OUT "none";
                }
                print OUT "\n";
            }
        }
        close $in;
        close OUT;
    }
}  

##########################################################################

if($options{tree_list}){
    my @iList;
    open(iLIST, "<", $options{tree_list}) or die "ERROR: Can not open dfam index list: $options{tree_list} because: $!\n";
    while(<iLIST>){
        chomp;
        push(@iList, $_);
    }
    close iLIST;
    my %chr_tree;
    foreach my $index (@iList){
        $index =~ /(\w+)\.tree$/;
        my $chr = $1;
        $chr_tree{$chr} = IntervalTree->new();
        $chr_tree{$chr}->treeFromFile($index);
    }
    foreach my $bam (@bam_list){
        open(my $in, "-|", "samtools view $bam") or die "Unable to open the input bam file: $bam because $!\n";    ## Read in input BAM
        $bam =~ m/(\w+?)\.bam$/;
        open(OUT, ">", "$options{output_dir}./$1\.dfam_overlap") or die "Unable to open output $options{output_dir}/$1.dfam_overlap because $!\n";
        while(<$in>){
            chomp;
            my @line=split(/\t/,$_);
            if($line[2] =~/chrM/){
               print OUT "$line[0]\tchrM\n";
               next;
            }  
            print OUT "$line[0]\t";
            my @overlaps = $chr_tree{$line[2]}->searchInterval($line[3], $line[3]+length($line[9]));      ## Find the intervals for each read
            if(@overlaps){                                                                                ## If output, print it
                foreach my $p (@overlaps){
                    print OUT "$p->[2]\t";
                }
            } elsif (!$chr_tree{$line[2]}->searchInterval($line[3], $line[3]+length($line[9]))){
                print OUT "none";
            }
            print OUT "\n";
        }
        close $in;
        close OUT;
    }
}

##############################################################################

## Process with one dfam file
if($options{dfam_file}){
    ## Build Tree for dfam TE hits in hg19
    ## $tree = IntervalTree->new(); 
    open(my $dfam_file, "<", "$dfam_file") or die "Unable to open the dfam_file $dfam_file because $!\n";
    my %chr_tree;
    while(<$dfam_file>){
        chomp;
        my @f=split(/\t/,$_);                                                             
        my @hit = split(/\./,$f[15]);
        my $chr;
        my %seen;
        if($hit[0] =~ /^GL(\d+)/){                                                       ## Next couple lines convert Dfam chr nomenclature to hg19 nomenclature
            $chr = "chrUn_gl$1";                  
        } elsif ($hit[0] =~ /^(\d+|X|Y)/){
            $chr = "chr$1";
        }
        if(!$chr_tree{$chr}){                                                             ## While reading through the dfam file, if we haven't seen the chr before, start a new tree for that chr
            $chr_tree{$chr} = IntervalTree->new();
            $seen{$chr}++;
        } 
        if($chr_tree{$chr}){
            if($f[9] =~/\+/){                                                             ## If on the positive strand
                $chr_tree{$chr}->addInterval($f[1],$f[10],$f[11]);                                  ##Add dfam acc.#, start and stop to tree
            }
            if($f[9] =~ /\-/){                                                            ## If on the negative strand
                $chr_tree{$chr}->addInterval($f[1],$f[11],$f[10]);                                  ##Add dfam acc.#, start and stop to tree
            }
        }
    }       
    close $dfam_file;
    for my $keys (keys %chr_tree){
        $chr_tree{$keys}->buildTree();
    }
    ## Finished building tree, now take each bam read and figure out if they overlap
    foreach my $bam (@bam_list){
        open(my $in, "-|", "samtools view $bam") or die "Unable to open the input bam file: $bam because $!\n";    ## Read in input BAM
        $bam =~ m/(\w+?)\.bam$/;
        open(OUT, ">", "$options{output_dir}./$1\.dfam_overlap") or die "Unable to open output $options{output_dir}/$1.dfam_overlap because $!\n";
        while(<$in>){
            chomp;
            my @line=split(/\t/,$_);
#            if($line[2] =~/chrM/){
#              print OUT "$line[0]\tchrM\n";
#              next;
#           }
            if(!$seen{$line[2]}){
               print OUT "$line[0]\t$line[2]\tUnable to search for this ...\n";
               next;
            }  
            my @overlaps = $chr_tree{$line[2]}->searchInterval($line[3], $line[3]+length($line[9]));      ## Find the intervals for each read
            print OUT "$line[0]\t";
            if(@overlaps){                                                                                ## If output, print it
                foreach my $p (@overlaps){
                    print OUT "$p->[2]\t";
                }
            } elsif (!$chr_tree{$line[2]}->searchInterval($line[3], $line[3]+length($line[9]))){
                print OUT "none";
            }
            print OUT "\n";
        }
        close $in;
        close OUT;
    }
}

##########################################################################

## Process with dfam's split by chr.
if($options{dfam_dir}){
    foreach my $bam (@bam_list){
        open(my $in, "-|", "samtools view $bam") or die "Unable to open the input bam file: $bam because $!\n";    ## Read in input BAM
        $bam =~ m/(\w+?)\.bam$/;
        open(OUT, ">", "$options{output_dir}./$1.dfam_overlap") or die "Unable to open output $options{output_dir}/$1.dfam_overlap because $!\n";
        while(<$in>){
            chomp;
            my @line=split(/\t/,$_);                                                           ## bam line split
            if($line[2] =~/chrM/){
               print OUT "$line[0]\tchrM\n";
               next;
            }  
            if(!$current_chr{$line[2]}){                                                       ##If I haven't seen this chr yet
                $tree = IntervalTree->new();                                                   ##clear tree
                open(my $dfam, "<", "$dfam_dir/$line[2].Dfam_hits") or next "Unable to open Dfam_hits file for $line[2]. Skipping\n";                           ##open dfam hits file to build tree for the new chr
                while(<$dfam>){ 
                    chomp;
                    my @f=split(/\t/,$_);
                    if($f[9] =~/\+/){                                                             ## If on the positive strand
                        $tree->addInterval($f[1],$f[10],$f[11]);                                  ##Add dfam acc.#, start and stop to tree
                    }
                    if($f[9] =~ /\-/){                                                            ## If on the negative strand
                        $tree->addInterval($f[1],$f[11],$f[10]);                                  ##Add dfam acc.#, start and stop to tree
                    }
                }
                close $dfam;
                $tree->buildTree();
                $current_chr{$line[2]}=1;                                                     ## 
            }
            if($current_chr{$line[2]}){                                                       ## If we already have the tree built for this chr 
               my @overlaps = $tree->searchInterval($line[3], $line[3]+length($line[9]));      ## Find the intervals for each read
                print OUT "$line[0]\t";
                if(@overlaps){                                                                ## If output, print it
                    foreach my $p (@overlaps){
                        print OUT "$p->[2]\t";
                    }
                } elsif (!$current_chr{$line[2]}->searchInterval($line[3], $line[3]+length($line[9]))){
                    print OUT "none";
                }
                print OUT "\n";
            } else {
                print STDERR "$line[2] doesn't have a dfam hits file.\n";
            }
        }
        close $in;
        close OUT;
    }
}

##########################################################################
   

##########################################################################
    
__END__
