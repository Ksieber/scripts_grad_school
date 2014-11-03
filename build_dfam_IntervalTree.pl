#!/usr/local/bin/perl
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use IntervalTree;

my %options;
my $results = GetOptions (\%options,
                          'input=s',
                          'input_list=s',
                          'header=s',
                          'output_dir=s',
                          'help',
);

if($options{help}){
    die "HELP: This script will build an IntervalTree (input) and write an ~index to disk (output).
\t--input=\t\tFile to build the interval tree on. Currently set up to parse Dfam_hits file. May need to change. 
\t--input_list=\t\tList of files to index.
\t--header=\t\tThis is either a bam file, or a samtools view header from a bam with all the possible chr lengths. 
\t--output_dir=s\t\tDirectory for the indexed filesto be writen.\n";
}

my @input_list;
my @header;
my %chr_max;

if(!$options{output_dir}){
    die "ERROR: Must use --output_dir=\n";
}

if(!$options{header}){
    die "Must give a header file or bam to get the chr lengths from. Use --header= (foo.header or foo.bam)\n";
}

if($options{input_list}){
    open(LIST, "<", "$options{input_list}") or die "ERROR: Can not open input_list: $options{input_list} because: $!\n";
    while(<LIST>){
        chomp;
        push(@input_list,$_);
    }
    close LIST or die "ERROR: Can not close input_list: $options{input_list} because: $!\n";
}

if($options{input}){
    push(@input_list, $options{input});
}



if($options{header}=~/\.bam$/){
    @header=`samtools view -H $options{header}`;
} else {
    open(HEADER, "<", $options{header}) or die "ERROR: Can not open $options{header} for reading because: $!\n";
    push(@header, $_);
}
foreach my $header (@header){
    $header =~ /\@SQ\s+SN:(chr\w+|X|Y|\_)\sLN:(\d+)/;                           ## Grab chr and max length for that chr
    $chr_max{$1}=$2;
}

foreach my $input (@input_list){
    $input =~ /(\w+)\.\w+$/;
    my $output = "$options{output_dir}/$1\.tree";
    my $tree = IntervalTree->new();
    #    my $max = 100;
    #    $tree = IntervalTree->max();
    #    $tree->"max" = $max;
#    $tree{"max"}= 100;
    open(IN, "<", "$input") or die "ERROR: Can not open input: $input because: $!\n";
    while(<IN>){
        chomp;
        my @f=split(/\t/,$_);                                                             
        my @hit = split(/\./,$f[15]);
        if($f[9] =~/\+/){                                                             ## If on the positive strand
            $tree->addInterval($f[1],$f[10],$f[11]);                                  ##Add dfam acc.#, start and stop to tree
        }
        if($f[9] =~ /\-/){                                                            ## If on the negative strand
            $tree->addInterval($f[1],$f[11],$f[10]);                                  ##Add dfam acc.#, start and stop to tree
        }
    }
    close IN or die "ERROR: Can not close $input: $input because: $!\n";
    $tree->buildTree();
    $tree->treeToFile($output);
}

__END__
 #my $chr;
#        if($hit[0] =~ /^GL(\d+)/){                                                    ## Next couple lines convert Dfam chr nomenclature to hg19 nomenclature
            $chr = "chrUn_gl$1";                  
        } elsif ($hit[0] =~ /^(\d+|X|Y)/){
            $chr = "chr$1";
#        }
