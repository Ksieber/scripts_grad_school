package setup_output;
use strict;
use File::Basename;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( setup_output );
## This script will setup output info
## Suggested use:
## my $out = &setup_output({input => $input});
## open(OUT,"> $out->{$file}->{dir}/foo.txt");
## open(OUT,"> $out->{$file}->{name}.txt");
## ARGS
## input => ref array of input files to process
## output_dir => directory for main output
## subdirs => <0|1> [0] 1=Makes a subdirectory in the output_dir foreach input


sub setup_output {
    my %output;
    my $outname;
    my $config = shift;
    if(!$config->{input}){die "Must pass &setup_output an input =>\n";}
    my $in = $config->{input};
    my $subdirs = $config->{subdirs} ? $config->{subdirs} : 0;
    my @suffix = qw( .bam .txt .srt.bam .list _\d+.fastq .fastq );
    foreach my $file (@$in){
        my ($fn,$path,$suf)=fileparse($file,@suffix);
        my $outdir;
        if($subdirs==1){
            $outdir = $config->{output_dir} ? "$config->{output_dir}$fn/" : "$path/$fn/";
        } else {
            $outdir = $config->{output_dir} ? $config->{output_dir} : $path; 
        }
        `mkdir -p $outdir`;
        $outname = "$outdir$fn";
        $output{$file}=({name => $outname, dir => $outdir, fn => $fn});
    }
    return \%output;
}


1;
