#!/usr/local/bin/perl -w
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;
my %options;
my $results = GetOptions (\%options,
        'sra_files=s',
        'output_dir=s',
        'decrypt_script=s',
        'threads=s',
        'rm',
        'help',
        );

if($options{help}){die
    "--sra_files=
        --output_dir=
        --decrypt_script=
        --threads
        --rm\n";
}

open IN, "<$options{sra_files}" or die "Couldn't open $options{sra_files}\n";
while(<IN>) {
    chomp;
    &process_file($_);
}
sub process_file {
    my $sra_file = shift;
    my @files;
    if($sra_file) {
        my($name,$path,$suff) = fileparse($sra_file,'.lite.sra.ncbi_enc');
        print `scp diagmaster.igs.umaryland.edu:$sra_file $options{output_dir}`;
        print STDERR "$options{decrypt_script} $options{output_dir}/$name$suff -remove-encrypted\n";
        print `$options{decrypt_script} $options{output_dir}/$name$suff -remove-encrypted`;
        my $file = "$options{output_dir}/$name.lite.sra";
        $file =~ s/\/\//\//g;
        `/usr/local/packages/sratoolkit.2.1.8/fastq-dump --split-3 -O $options{output_dir} $file`;
        @files = `find $options{output_dir} -name '$name*.fastq'`;
    }

    map{chomp;}@files;
    my @srted = sort @files;
#print STDERR join("\t",@srted);
    my $fastq1 = $srted[0];
    my $fastq2 = $srted[1];
    if(scalar @srted <2) {
        die "Only had 1 fastq file\n";
    }

    $fastq1 =~ /(\w+)_1|2\.fastq/;
    my $pref=$1;
    my $t= $options{threads} ? "$options{threads}" : "1";
    `perl /home/ksieber/scripts/BWA_aligner.pl --fastq_1=$fastq1 --fastq_2=$fastq2 --ref=/local/projects-t3/HLGT/references/hg19/rna.fa --output_prefix=$pref --output_dir=$options{output_dir} --t=$t --bam_output --sort_index_bams --cleanup --mapped_only 2>$options{output_dir}/$pref\_error.log`;
    if ($? != 0 ){print STDERR "ERROR:$. Unable to align the input correctly?\n";}
    `samtools view $options{output_dir}/$pref.srt.bam 'gi|98986444|ref|NM_004363.2|' -bo $options{output_dir}/$pref\_ceacam5_rna_region.bam`;
    `samtools view $options{output_dir}/$pref.srt.bam 'gi|72255578|ref|NM_021103.3|' -bo $options{output_dir}/$pref\_tmsb10_rna_region.bam`;
    `samtools view $options{output_dir}/$pref.srt.bam 'gi|168480144|ref|NM_001101.3|' -bo $options{output_dir}/$pref\_actb_rna_region.bam`;
    if ($? != 0 ){print STDERR "ERROR:$. Unable to get reads from the regions of intrest. Either no reads to pull or error...\n";}
    if($options{rm} && $? == 0){
        `rm $options{output_dir}/$pref\.srt.ba*`;
        `rm $options{output_dir}/$pref\_lrg_insert.histogram`;
        `rm $options{output_dir}/$pref\_std_insert.histogram`;
        `rm $options{output_dir}/$pref\_1.fastq`;
        `rm $options{output_dir}/$pref\_2.fastq`;
        `rm $options{output_dir}/$pref\.lite.sra`;
    }
    print STDERR "Finished processing $sra_file\n";
}
