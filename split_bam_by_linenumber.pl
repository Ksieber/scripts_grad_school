#!/usr/bin/perl

=head1 NAME

lgtseek_split_bam.pl

=head1 SYNOPSIS

Splits input bams into chunks of a certain size.

=head1 DESCRIPTION

Splits input bams into chunks of a certain size.

=head1 AUTHOR - David R. Riley

e-mail: driley@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut
use lib '/home/ksieber/scripts/lgtseek-master/lib/';
use strict;
use Scalar::Util qw(reftype);
use Qsub;
use LGTSeek;
use File::Basename;
use setup_input;
use setup_output;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
my $results = GetOptions (\%options,
                          'input_list=s',
                          'input=s',
                          'line_count=s',
                          'Qsub',
                          'output_dir=s',
                          'subdirs',
                          'samtools_bin=s',
                          'ergatis_dir=s',
                          'output_list=s',
                          'bin_dir=s',
                          'help|h'
                          );

if($options{help}){
   die "Help: This script will split bams into smaller bams based on the number of lines per file. 
      --input=         <BAM>
      --input_list=    <LIST>
      --line_count=    <lines per bam> [250,000,000]
      --Qsub           Qsub the split command. Great for processing many bams at once.
      --output_dir=    Directory for output. 
      --subdirs        Make a directory within the output_dir to place the output. 
      --output_list=   Name the list of output files created.
      --bin_dir=       Directory where LGTSeek.pm is stored.
      --help\n";
}

# Take care of the inputs
if(!$options{input} && !$options{input_list}){die "Error: Must give input with --input=<BAM> or --input_list=<LIST of bams>\n";}
my $bin_dir = $options{bin_dir} ? $options{bin_dir} : '/local/projects-t3/HLGT/scripts/lgtseek-master/lib/';

my $ergatis_dir = $options{ergatis_dir} ? $options{ergatis_dir} :'/local/projects/ergatis/package-driley/bin/';

my $samtools_bin = $options{samtools_bin} ? $options{samtools_bin} : 'samtools';

if(!$options{line_count}){$options{line_count}=250000000};

# Create an lgtseek object
my $lgtseek = LGTSeek->new({
    bin_dir => $bin_dir,
    output_dir => $options{output_dir},
    ergatis_bin => $ergatis_dir,
    samtools_bin => $samtools_bin,
    paired_end => 1
});

## New Method: Handle input and output_dirs
my $input = setup_input();
my $out_dir = setup_output($input);


# Open a list file to write the output bams to
my $olistfh;
if($options{output_list}) {
    open($olistfh, ">$options{output_list}") or die "Unable to open $options{output_list}\n";
}

# Loop over the bam files and split them
foreach my $file (@$input) {
    if($options{Qsub}){
       my ($fn,$path,$sufx)=fileparse($file,".bam");
       if($path=~/^\.\/$/){ 
          die "Error: Use full path names for Qsub.\n";
       }
       open(SH,">","$out_dir->{$file}/$fn\_split.sh") || die "Error: Can't open shell scipt: $out_dir->{$file}/$fn\_split.sh because: $!\n";
       print SH "perl ~/scripts/split_bam_by_linenumber.pl --input=$file --line_count=$options{line_count} --output_dir=$out_dir->{$file}";
       close SH;
       Qsub("$out_dir->{$file}/$fn\_split.sh");
       next;
    }
    my $split_files = $lgtseek->splitBam({
        input => $file,
        seqs_per_file => $options{line_count},
        output_dir => $out_dir->{$file},
    });
    
    # Write out the filenames to the list file.
    if($olistfh) {
        map {
            print $olistfh "$_\n";
        }@$split_files
    }
}


__END__
