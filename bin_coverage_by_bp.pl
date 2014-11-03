#!/usr/local/bin/perl -w
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
      'vcf=s',
      'output_prefix=s',
      'output_dir=s',
      'help|h',
      );

if($options{help}){
   die "Help: This script will count the number of base positions with X coverage.
      \t--vcf=\t\tInput mpileup file. REQUIRED.
      \t--output_prefix=\tPrefix for the output. Optional.
      \t--output_dir=t\tDirectory for the output. Optional.\n";
}

if(!$options{vcf}){die "ERROR: Must give an input file. Use --vcf= name of the input mpileup file.\n";}
my %cov;
$options{vcf} =~ m/(\w+)\.\w+$/;
my $out = defined $options{output_prefix} ? "$options{output_prefix}" : "$1";
open(OUT, ">", "$options{output_dir}./$out\.bin_cov") or die "ERROR: Unable to open output: $options{output_dir}/$out\.bin_cov because: $!\n";

open(IN, "<", "$options{vcf}") or die "ERROR: Unable to open input vcf: $options{vcf} because: $!\n";
while(<IN>){
   my @f=split(/\t/);
   if($f[3] eq 1){$cov{1}++;}
   if($f[3] eq 2){$cov{2}++;}
   if($f[3] eq 3){$cov{3}++;}
   if($f[3] eq 4){$cov{4}++;}
   if($f[3] eq 5){$cov{5}++;}
   if($f[3] eq 6){$cov{6}++;}
   if($f[3] eq 7){$cov{7}++;}
   if($f[3] eq 8){$cov{8}++;}
   if($f[3] eq 9){$cov{10}++;}
   if($f[3] > 9 && $f[3] < 20){$cov{10}++;}
   if($f[3] > 19 && $f[3] < 50){$cov{20}++;}
   if($f[3] > 49 && $f[3] < 100){$cov{50}++;}
   if($f[3] > 99){$cov{100}++;}
}
close IN or die "ERROR: Unable to close input vcf: $options{vcf} because: $!\n";

print OUT "Coverage\tNumber of bp\n";
foreach my $key (sort keys %cov){
   print OUT "$key\t$cov{$key}\n"
}
close OUT or die "ERROR: Unable to close output: $options{output_dir}/$out\.bin_cov\n";

__END__
