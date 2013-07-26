#!/usr/bin/perl 
use warnings;
use strict; 
use Bio::SearchIO;
use Data::Dumper;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options = ();
my $results = GetOptions (\%options,
              'input=s',
              'output_dir=s',
              'output_name=s',
              'input_format',
              'max_eval=s',
              'min_length=s',
              'help|h');

if($options{help}){ die 
   "Help: This script will take a blast report and return each hits' BEST HSP in fasta format.
      --input=             Input blast report. 
      --input_format=      Format of the blast report. Options: [BLAST], PSIBLAST, PSITBLASTN, RPSBLAST, WUBLAST, bl2seq, WU-BLAST, BLASTZ, BLAT
      --output_dir=        Directory for output. Defaults to cwd. 
      --output_name=       Name of the output fasta. Defaults to \"input\"\_besthits.fa.
      --max_eval=          Max Evalue returned. Defaults to 10.
      --min_length=        Min HSP total length. Defaults to 20.\n";
}

######## Check Input ########
my $input = $options{input} ? $options{input} : $ARGV[0];
if(!$input){die "Must give an input with --input=<FILE> or with ARGV[0].\n";}
my ($fn,$path,$suf)=fileparse($input,".txt");
$path=~s/\/$//;
my $input_format = $options{input_format} ? $options{input_format} : 'blast';
my $max_eval = $options{max_eval} ? $options{max_eval} : '10';
my $min_length = $options{min_length} ? $options{min_length} : '20';

######### Set output #########
my $out_dir = $options{output_dir} ? $options{output_dir} : $path;
my $out_file = $options{output_prefix} ? $options{output_prefix} : "$fn\_besthits.fa";
my $out = "$out_dir/$out_file";
open(OUT,">","$out") || die "Can't open output: $out because: $!\n";

my $search = Bio::SearchIO->new(-format => $input_format,
								-file => $input,
								-best => $options{best},
								-signif => $max_eval,
);
									
while( my $result = $search->next_result ){
	while (my $hit = $result->next_hit ){
		while( my $hsp = $hit->next_hsp ){
			if($hsp->evalue <= $max_eval ){
            if($hsp->length('total')>=$min_length){
				   if($hsp->rank == 1){
                  my $seq=$hsp->hit_string;
                  $seq =~ s/-//g;
					   print OUT "\>".$hit->name." ".$hit->description."\n$seq\n";
               }
				}
			}
		}
	}
}
		
	
