#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use run_cmd;
use setup_input;
use File::HomeDir;
use Carp;
use mk_dir;
$Carp::MaxArgLen = 0;
use XML::Simple qw(:strict);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'input_list|I=s', 'output_dir=s', 'sub_dirs=i', 'help|?', 'Qsub|q=i', 'Qsub_iterate=i', 'sub_mem=s', 'sub_sleep=s', )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) { &help; }
if ( !$options{input} && !$options{input_list} ) { die "Error: Must pass the full path of a bam(s) to --input or --input_list\n"; }

$options{sub_mem}   = $options{sub_mem}   ? $options{sub_mem}   : "100M";
$options{sub_sleep} = $options{sub_sleep} ? $options{sub_sleep} : ".1";

if ( $options{Qsub} or $options{Qsub_iterate} ) { Qsub_script( \%options ); }

my $inputs = setup_input( \%options );

foreach my $input (@$inputs) {
    my ( $fn, $path, $suff ) = fileparse( $input, "\.bam" );
    my @split_path  = split( /\//, $path );
    my $analysis_id = $split_path[-1];
    my $output_dir  = $options{output_dir} ? $options{output_dir} : File::HomeDir->my_home;
    if ( $options{sub_dirs} == 1 ) { $output_dir = $output_dir . "$analysis_id"; }
    mk_dir($output_dir);
    my $xml = "$output_dir/cgquery.xml";

    # First, pull cgquery xml file from cghub
    my $cmd              = "/home/ksieber/lib/Python-2.7.7/bin/python2.7 /home/ksieber/bin/cgquery \"analysis_id\=$analysis_id\" -o $xml";
    my $cgquery_retry    = 1;
    my $cgquery_attempts = 0;
    my $retry_pause      = 3;
    while ( $cgquery_retry == 1 && $cgquery_attempts <= 3 ) {
        my $cgquery_exec_return = `$cmd`;
        if ($?) {
            print STDERR "***Error*** :: $cmd :: failed with message: $cgquery_exec_return :: $?. Will now retry cgquery.\n";
            $cgquery_attempts += 1;
            $retry_pause = $retry_pause * 5;
            sleep $retry_pause;
        }
        else {
            $cgquery_retry = 0;
        }
    }

    # Second, parse md5 value for the bam
    my $cghub_md5;
    my $ref = XMLin( $xml, KeyAttr => { Result => 'id' }, ForceArray => [ 'Result', 'file' ], GroupTags => { files => 'file' } )
        or confess "&XMLin is unable to read in the cgquery.xml file: $xml because: $!\n";
    foreach my $sample ( sort keys %{ $ref->{'Result'} } ) {
        foreach my $file ( @{ $ref->{'Result'}->{$sample}->{files} } ) {
            if ( $file->{filename} =~ /\.bam$/ ) {
                $cghub_md5 = $file->{checksum}->{content};
            }
        }
    }
    run_cmd("rm $output_dir/cgquery.xml");

    # Third, calculate md5
    my $calc_md5 = run_cmd("md5sum -b $input | cut -f1 -d \" \"");
    chomp($calc_md5);

    # Lastly, compare the calculated md5 to the cghub md5
    if ( $calc_md5 eq $cghub_md5 ) {
        print STDOUT "Pass: $input\n";
        next;
    }
    elsif ( $calc_md5 ne $cghub_md5 ) {
        print STDOUT "*** WARNING *** md5 doesn't match cghub for: $input | Calculated md5: $calc_md5 | CGQuery md5: $cghub_md5\n";
    }
}

sub help {
    die "\nThis script will check the md5 checksum for a TCGA bam.
    --input=        /full/path/to/tcga_analysis_id/file.bam
    --input_list=   /list/of/bams
    --output_dir=       [~]
    --Qsub=
    --Qsub_iterate=\n\n";
}
