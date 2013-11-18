#!/usr/bin/perl
=head1  NAME 

lgt_overlap_finder.pl

=head1  SYNOPSIS

    lgt_overlap_finder.pl 
        [--genbank_file=/a/file.gbk || --genbank_dir=/a/dir/of/gbks/]
        --input=/kellys/crappy/file.bam
        --output_dir=/kellys/face/
        [-samtools_path=/path/to/samtools]

=head1   DESCRIPTION


=head1 INPUT



=head1 OUTPUT


=head1 CONTACT

    David Riley
    driley@som.umaryland.edu

=cut

use strict;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Basename;
use Data::Dumper;
use Bio::SeqIO;
use IntervalTree;
use Bio::DB::EUtilities;
use File::Find;

$|++;

my %options = ();
my $results = GetOptions (
    \%options,
    'genbank_file=s',
    'genbank_dir=s',
    'input=s',
    'output_dir=s',
    'samtools_path:s',
    'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    die "This script will find the genes that lgt reads overlap with.
    --input=            <BAM>
    --genbank_file=     
    --genbank_dir=  
    --output_dir=
    --help\n";
}
if(!$options{input}) {
    pod2usage( {-exitval=>0, -verbose => 1, -output => \*STDOUT} );

}
my $samtoolsbin = $options{samtools_path} ? "$options{samtools_path}/samtools" : "/usr/local/bin/samtools";
# If a genbank file hasn't been passed in then we'll figure out the accession
# looking at the bam header.
if($options{genbank_dir}) {
    &find_gb_file();
}
elsif(!$options{genbank_file}) {
    &get_gb_file();
}


my $tree = IntervalTree->new();
#First read in the genbank file
print STDERR "Reading genbank file\n";
&read_genbank();
print STDERR "Building tree\n";
$tree->buildTree();

print STDERR "Processing bam\n";
&read_bam();

sub get_ref_id {
    my @lines = `$samtoolsbin view -H $options{input}`;
    if($?) {
        die "Had a problem $?\n";
    }
    my $id;
    map {
        if(/^\@SQ\s+SN:gi\|\d+\|\w+\|(\S+)\|/) {
            $id = $1;
        }
    }@lines;
    if(!$id) {
        die "Unable to pull the id from $options{input}\n";
    }
    return $id;
}

sub find_gb_file {
    my $file;

    my $id = &get_ref_id();
    $id =~ s/\..*//;
    print STDERR "Looking for $id.gbk in $options{genbank_dir}\n";
    find(sub { 
        if($File::Find::name =~ /$id.gbk/) { 
            $file = $File::Find::name;
       }},$options{genbank_dir});
    $options{genbank_file} = $file;
}

sub get_gb_file {

    my $id = &get_ref_id();

    my $efetch = Bio::DB::EUtilities->new(
                                          -db => 'nucleotide',
                                          -id => $id,
                                          -rettype => 'gb',
                                          -retmode => 'text',
                                          );
    
    open OUT, ">$options{output_dir}/$id.gbk" or die "Unable to open $options{output_dir}/$id.gbk\n";
    print OUT $efetch->get_Response->content;
    close OUT;

    $options{genbank_file} = "$options{output_dir}/$id.gbk";
}

sub read_bam {
    open(my $handle,"-|", "$samtoolsbin view $options{input}");

    while(my $line = <$handle>) {
        my @fields = split(/\t/,$line);
        my $flag = &parseFlag($fields[1]);
        if($flag->{'query_mapped'}) {
            my @overlaps = $tree->searchInterval($fields[3], $fields[3]+$fields[9]);
            if(@overlaps) {
                foreach my $p (@overlaps) {print "$fields[0]\t$p->[2]->{primary}\t$p->[2]->{product}\t$p->[2]->{locus}\n";}
            } else {
                print "No overlaps found for $fields[0]\n";
            }
        }
    }
}

sub read_genbank {

    my $file = "$options{genbank_file}";

    my $stream = Bio::SeqIO->new(-file => $file,
                              -format => 'GenBank');
    
    while(my $seq = $stream->next_seq()) {
        my @feats = $seq->get_all_SeqFeatures();

        foreach my $feat (@feats) {
            my $primary = $feat->primary_tag();
            # Basically we're only looking for things with a product.
            # This includes rRNA,tRNA and CDS I believe
            if($feat->has_tag('product') && $feat->has_tag('locus_tag')) {
                my @vals = $feat->get_tag_values('product');
                my @locus = $feat->get_tag_values('locus_tag');
                my $start = $feat->start;
                my $stop = $feat->end;
                my $obj = {
                    'primary' => $primary,
                    'product' => $vals[0],
                    'locus' => $locus[0],
                    'start' => $feat->start,
                    'end' => $feat->end
                };
                $tree->addInterval($obj,$feat->start,$feat->end);
            }
        }  
    }
}

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}

sub parseFlag {
    my $flag = shift;
    my $rawbin = dec2bin($flag);
    my $rev = scalar $rawbin;
    if($rev eq $rawbin) {
        #    print "ERROR $rev $rawbin\n";
    }
    my $bin = sprintf("%011d", $rev);
    my $final_bin = reverse $bin;
    my $prop = substr($final_bin, 1, 1);
    my $qmap = substr($final_bin, 2, 1);
    my $mmap = substr($final_bin, 3, 1);
    my $qstrand = substr($final_bin, 4, 1);
    my $mstrand = substr($final_bin, 5, 1);

    return {
        'query_mapped' => !$qmap,
        'mate_mapped' => !$mmap
    };
}
