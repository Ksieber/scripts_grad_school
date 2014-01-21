#!/usr/local/bin/perl 

use warnings;
use strict;
use run_cmd;
use Carp;
$Carp::MaxArgLen = 0;
use File::Basename;
use lib qw(/local/projects-t3/HLGT/scripts/lgtseek/lib/);
use LGTSeek;

if ($#ARGV != 2) {
   print "\nNeed 3 parameters: 
            Input: BAM or mpileup
            window size
	    	distance between windows\n";
   exit;
}
my $lgtseq = LGTSeek->new2;
my ($file,$window,$distance) = @ARGV;
my %hash;

#my $A = defined $options{A} ? $options{A} : "1";
#my $d = defined $options{d} ? $options{d} : 100000;
#my $region = defined $options{region} ? "'$options{region}'" : undef;
#my $sort = defined $options{sort} ? "$options{sort}" : "0";
my $region = undef;
my $sort = "0";
my $A = "1";
my $d = "100000";
my $view = "-hu";              ## Default
my $mpileup = "-Ad $d";    	   ## Default
if($A==0){ $mpileup = "-d $d"; }

my $fh = &open_input($file);
while(<$fh>){
    chomp;
    my ($chr,$coord,$value)=(split)[0,1,3];
    $hash{$chr}{$coord}=$value;
}
close ($fh) || die "Cannot close $file : $!\n";

foreach my $chr (keys %hash) {
#	print STDERR "$chr\n";
	foreach my $pos (keys %{$hash{$chr}}){
#		print STDERR "$pos\n";
    	my $a++;
    	my $sum=0;
    	my $e=1;

    	for (my $i=1; $i<$window; $i++) {
			if (exists $hash{$chr}{$i}) {
	    		$sum=$sum+$hash{$chr}{$pos+$i};
	    		$e++;
	    	}
    	}
    	my $average=$sum/$e;
    	my $test=$distance;
    	if ($a>$test) {
        	print "$chr\t$pos\t$average\n";
			$test=$test+$distance;
    	}
	}
}

sub open_input{
	my $raw_input = shift;
	my ($fn,$path,$suffix)=fileparse($raw_input,(@{$lgtseq->{mpileup_suffix_list}},@{$lgtseq->{bam_suffix_list}}));
	
	if($lgtseq->empty_chk({input => $raw_input})==1){confess "ERROR: Input file:$raw_input is empty.";}
	## Filehandle to return.
	my $fh;
	## Mpileup input
	map {
		if($suffix eq $_){
			open($fh, "<","$raw_input") || confess "ERROR: Can't open input: $raw_input because: $!\n"; 
			return $fh;
		}
	} @{$lgtseq->{mpileup_suffix_list}};

	## Position sorted bam
	my @psort_bams = ('_pos-sort.bam','_psort.bam','srt.bam');
	map {
		if($suffix eq $_){
			open($fh, "-|","samtools view $view $raw_input $region | samtools mpileup $mpileup - ") || confess "ERROR: Can't open input: $raw_input because: $!\n"; 
			return $fh;
		}
	} @psort_bams;


	## Name sorted bam || --sort=1
	map {
		if($suffix eq $_ || $sort==1){
			print STDERR "Sorting input bam . . .\n";
			open($fh, "-|","samtools sort -m 5G -o $raw_input - | samtools view $view - $region | samtools mpileup $mpileup - ") || confess "ERROR: Can't open input: $raw_input because: $!\n"; 
			return $fh;
		}
	}@{$lgtseq->{bam_suffix_list}};
}

__END__
