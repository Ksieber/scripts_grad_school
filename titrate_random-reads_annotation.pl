#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/
use strict;
use warnings;
use lib ( '/home/ksieber/scripts/', '/local/projects-t3/HLGT/scripts/lgtseek/lib' );
use run_cmd;
use Cwd;
use mk_dir;
use bwa;

if ( !@ARGV ) { &help; }

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results
    = GetOptions( \%options, 'ref|r=s', 'annotation|a=s', 'num|n=i', 'read_length|l=i', 'SE', 'output_dir|o=s', 'output_prefix|p=s', 'bwa=s', 'wgsim=s', 'mem=i', 'read_distribution_script|rds=s',
    'threads|t=i', 'help|?' )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( !defined $options{ref} )        { die "Error: Must provide: --ref\n"; }
if ( !defined $options{annotation} ) { die "Error: Must provide: --annotation}\n"; }
if ( !-e $options{ref} )             { die "Error: The reference provided doesn't exit: $options{ref}\n"; }
if ( !-e $options{annotation} ) { die "Error: The annotation file provided doesn't exist: $options{annotation}\n"; }

my $bwa                      = defined $options{bwa}                      ? $options{bwa}                      : "/home/ksieber/bin/bwa";
my $wgsim                    = defined $options{wgsim}                    ? $options{wgsim}                    : "/home/ksieber/bin/wgsim";
my $read_distribution_script = defined $options{read_distribution_script} ? $options{read_distribution_script} : "/home/ksieber/bin/read_distribution.py";
my $read_length              = defined $options{read_length}              ? $options{read_length}              : "100";
my $iterations               = defined $options{num}                      ? $options{num}                      : "1000";
my $pe                       = defined $options{SE}                       ? "0"                                : 1;
my $output_dir               = defined $options{output_dir}               ? $options{output_dir}               : getcwd;
my $output_prefix            = defined $options{output_prefix}            ? $options{output_prefix}            : "titrate_random_annotation";
my $threads                  = defined $options{threads}                  ? $options{threads}                  : 8;
my $bwa_threads      = ( ( $threads - 2 ) > 2 )            ? ( $threads - 2 )            : 2;
my $samtools_threads = ( ( $threads - $bwa_threads ) > 2 ) ? ( $threads - $bwa_threads ) : 2;
my $out              = $output_dir . $output_prefix;
my $ref              = $options{ref};
my $annotation       = $options{annotation};

mk_dir($output_dir);
chdir($output_dir);

my $data;

for ( my $n = 1; $n <= $iterations; $n++ ) {
    my $tmp_dir = "$out/tmp_dir_$n";
    mk_dir($tmp_dir);

    if ( $pe == 1 ) {
        run_cmd("$wgsim -1 $read_length -2 $read_length $ref $tmp_dir/random_1.fq $tmp_dir/random_2.fq");
        bwa_mem(
            "$tmp_dir/random_1.fq",
            "$tmp_dir/random_2.fq",
            $ref,
            {   output_dir    => $tmp_dir,
                output_prefix => "test_$n",
                threads       => $bwa_threads
            }
        );
    }
    else {
        run_cmd("$wgsim -1 $read_length -2 $read_length $ref $tmp_dir/random_1.fq /dev/null");
        run_cmd("bwa mem -t $threads $ref $tmp_dir/random_1.fq | samtools view -@ $samtools_threads -S - -bo $tmp_dir/test_$n\.bam");
    }

    run_cmd("perl /home/ksieber/scripts/bam2srtbai.pl $tmp_dir/test_$n\.bam");
    run_cmd("$read_distribution_script -i $tmp_dir/test_$n\.psort.bam -r $annotation > $tmp_dir/overlap_$n.txt");

    open( IN, "<", "$tmp_dir/overlap_$n.txt" ) or die "Error: Couldn't open: $tmp_dir/overlap_$n.txt\n";
    while (<IN>) {
        next if ( $_ =~ /^Total/ );
        next if ( $_ =~ /^=/ );
        next if ( $_ =~ /^Group/ );
        my @f = split( /\s+/, $_ );
        push( @{ $data->{ $f[0] } }, $f[2] );
    }
    close IN;
    run_cmd("rm -rf $tmp_dir");
}

my $medians;
foreach my $key ( keys %{$data} ) {
    my $median = Math::NumberCruncher::Median( $data->{$key} );
    $medians->{$key} = $median;
}

foreach my $key ( sort { $medians->{$a} <=> $medians->{$b} } keys %{$medians} ) {
    print "$key\t$medians->{$key}\n";
}

sub help {
	die "This script isn't complete.\n";
}

