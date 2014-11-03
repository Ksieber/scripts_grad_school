#!/usr/local/bin/perl 

## Makes a sliding window calculation of a two column text file
## First column: position
## Second column: value to be averaged

if ($#ARGV != 2) {
   print "\nNeed 3 parameters: 
            file.txt (Space delimited file of values to be averaged)
            window size
	    distance between windows\n";
   exit;
}

($file,$window,$distance) = @ARGV;

open (FILE, "$file") || die "Cannot open $file : $!\n";
while (<FILE>) {
    chomp;
    ($coord,$value)=(split)[0,1];
    $hash{$coord}=$value;
}
close (FILE) || die "Cannot close $file : $!\n";

foreach $key (keys%hash) {
    $a++;
    $sum=0;
    $e=0;

    for ($i=0; $i<$window; $i++) {
	if (exists $hash{$i}) {
	    $sum=$sum+$hash{$key+$i};
	    $e++;
	    }
    }
    $average=$sum/$e;
    $test=$distance;
    if ($a>$test) {
        print "$key\t$average\n";
		$test=$test+$distance;
    }
}

__END__
