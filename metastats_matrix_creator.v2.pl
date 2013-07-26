#!/usr/local/bin/perl
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options;
my $results = GetOptions (\%options,
   	    'LGT_file=s',
	    'MM_file=s',
        'liberal_lca',
	    'help|h',
	    );

if ($options{help}){die "\nHELP: This script will take two files and create the appropriate input for metastats.  Output goes to STDOUT.
\n--LGT_file= the file with LGT read names and bacterial taxa
--MM_file= file with M_M read names and bacterial taxa
--liberal_lca   Analyze the lca for the liberal_lca of the second column in sam2lca.out

Note:Input= uniq.readname     taxa\n";}

##Variables
my %lgt_info;
my $lgt_subject;
my $lgt_bacteria;
my %mpd_info;
my $mpd_subject;
my $mpd_bacteria;
my %bacteria_unique;
my $NA=1;

#Make Hash of Hashes (hash{Subject}->Bacteria->Bacterial_count) for LGT data
if($options{LGT_file}){
    open (LGT, "<", $options{LGT_file}) || die "ERROR: Couldn't open $options{LGT_file} because: $!\n";
    while (<LGT>){
	   chomp;
	   my ($lgt_read,$lgt_bacteria)=split(/\t/, $_);
	   if ($lgt_read =~ m/^(\w+?)\./) {
	      $lgt_subject = $1;
      }
  	   $lgt_info{$lgt_subject}->{$lgt_bacteria}++;
    }
    close LGT || die "ERROR: Couldn't close $options{LGT_file} because: $!\n";
}


#Makes the same Hash of Hashes but for M_M data
open (MPD, "<", $options{MM_file}) || die "ERROR: Couldn't open $options{MM_file} because: $!\n";
while (<MPD>){
    chomp;
    if($options{liberal_lca}){
        my @F = split(/\t/, $_);
        $F[0] =~/^(\w+?)\./;   ##Grab SRR #
        $mpd_info{$1}->{$F[2]}++;
    } else {
        my ($mpd_read,$mpd_bacteria)=split(/\t/, $_);
        if ($mpd_read =~ m/^(\w+?)\./){
            $mpd_subject = $1;
        }
        $mpd_info{$mpd_subject}->{$mpd_bacteria}++;
    }
}
close MPD || die "ERROR: Couldn't close $options{MM_file} because: $!\n";


#Start creating metastats matrix.  
#First step is to create subject headers (twice, 1 for LGT 1 for MM)
#Also create "subject.list.tmp" with single subject/line to wc-l for count of # of subjects
#open (SUBJ, ">", "subject.list.tmp") || die "ERROR: Couldn't open subject.list.tmp because: $!\n";
if($options{LGT_file}){
    for my $k1 ( sort keys %lgt_info ) {
        if($k1 =~ /^\w+?/){
            print STDOUT "\tLGT$k1";
            #print SUBJ "$k1\n";
        } else {
            print STDOUT "\tFoo$NA";
            $NA++;
        }  
    }
    for my $k1 ( sort keys %lgt_info ) {
        if($k1 =~ /^\w+?/){
            print STDOUT "\tBacMM$k1";
        } else {
            print STDOUT "\tFoo$NA";
            $NA++;
        }
    }
} elsif ($options{MM_file}){
    for my $k1 ( sort keys %mpd_info ) {
        if($k1 =~ /^\w+?/){
            print STDOUT "\tBacMM$k1";
        } else {
            print STDOUT "\tFoo$NA";
            $NA++;
        }
    }
}

print STDOUT "\n";


#Create unique "list" of found bacteria in LGT 
if($options{LGT_file}){
    for my $k1B ( sort keys %lgt_info ) {
	   for my $k2 ( sort keys %{%lgt_info->{$k1B}} ) {
	      $bacteria_unique{$k2}++;
	   }
    }
} elsif ($options{MM_file}){
    for my $k1B ( sort keys %mpd_info ) {
	   for my $k2 ( sort keys %{%mpd_info->{$k1B}} ) {
	      $bacteria_unique{$k2}++;
	   }  
    }
}


#For each bacteria in LGT, print number of times it was seen in LGT subject first, then MM
for my $each_unique_bacteria (sort keys %bacteria_unique) {
    if ($each_unique_bacteria =~ /^(\w+?)/){
       print "$each_unique_bacteria\t";
    } else {
       print "Unknown\t";
	 }
    if($options{LGT_file}){
	   for my $k1C ( sort keys %lgt_info) {
	      if ($lgt_info{$k1C}->{$each_unique_bacteria}){
		      print STDOUT "$lgt_info{$k1C}->{$each_unique_bacteria}\t";
	      } else {print "0\t";}
	   }  
      for my $kD ( sort keys %lgt_info) {
	      if ($mpd_info{$kD}->{$each_unique_bacteria}){
		      print STDOUT "$mpd_info{$kD}->{$each_unique_bacteria}\t";
	      } else {print "0\t";}
	   }
	   print STDOUT "\n";
     }elsif($options{MM_file}){
	      for my $kD ( sort keys %mpd_info) {
	         if ($mpd_info{$kD}->{$each_unique_bacteria}){
		         print STDOUT "$mpd_info{$kD}->{$each_unique_bacteria}\t";
	         } else {print "0\t";}
	      }
	      print STDOUT "\n";	
    }
}

##`rm subject.list.tmp`;
