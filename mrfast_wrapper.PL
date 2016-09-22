#!/usr/bin/perl
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use strict;
use warnings;
use mk_dir;
use print_call;
use run_cmd;
use File::Basename;
if ( !@ARGV ) { &help; }
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options,       'fastq1=s',           'fastq2=s',     'fastq=s',     'min=i',           'max=i',             'output_dir|o=s', 'cleanup_none_VH_files=i',
    'auto_append=i', 'append_lane_number', 'append_lib=s', 'append_rg=s', 'append_sample=s', 'output_prefix|p=s', 'ref|r=s',        'subdirs=i',
    'sub_name=s',    'sub_mem=s',          'gz=i',         'help|?',      'Qsub|q=i',
) or die "Unrecognized command line option. Please try agian.\n";

if ( $options{help} )         { &help; }
if ( !defined $options{ref} ) { die "Error: Must use --ref=/path/to/ref.fa\n"; }
if ( !defined $options{min} or !defined $options{max} ) { die "Error: Must use --min=# and --max=# to insert size to use as discordant cutoff.\n"; }

my $interleaved = defined $options{fastq} ? 1 : 0;
my $clean_non_VH_files = defined $options{cleanup_none_VH_files} ? $options{cleanup_none_VH_files} : "1";

my $tmp_input = defined $options{fastq} ? $options{fastq} : $options{fastq1};
if ( !-e $tmp_input ) { die "Error: The input file doesn't exist: $tmp_input\n"; }
my ( $tmp_fn, $tmp_dir, $tmp_suf ) = fileparse( $tmp_input, qr/(_R[12]{1})?(\_\d+)?\.f\w{0,3}q(.gz)?/ );
my $out_dir    = defined $options{output_dir}    ? $options{output_dir}    : $tmp_dir;
my $out_prefix = defined $options{output_prefix} ? $options{output_prefix} : $tmp_fn;
if ( defined $options{subdirs} and $options{subdirs} == 1 ) { $out_dir = "$out_dir$out_prefix/"; }
mk_dir($out_dir);
$tmp_suf =~ /(_R([12]){1})?(\_\d+)?(\.f\w{0,3}q)(.gz)?/;
my ( $read_name, $read_number, $split_number, $fastq_suffix_name, $gz_suf ) = ( $1, $2, $3, $4, $5 );
my $lane_number;
if ( $tmp_fn =~ /\_(L\d{3})$/ ) { $lane_number = $1; }
my $file_is_gzipped         = ( $gz_suf eq '.gz' ) ? "1"          : "0";
my $override_gzipped_status = defined $options{gz} ? $options{gz} : $file_is_gzipped;
my $output                  = "$out_dir$out_prefix$split_number\.sam";

if ( defined $options{Qsub} and $options{Qsub} == 1 ) {
    $options{output_dir} = $out_dir;
    $options{subdirs}    = 0;
    $options{sub_name}   = defined $options{sub_name} ? $options{sub_name} : "mrfastALN";
    $options{sub_mem}    = defined $options{sub_mem} ? $options{sub_mem} : "5G";
    Qsub_script( \%options );

}

print_call( \%options );

my $mrfast_cmd = "/home/ksieber/bin/mrfast --search $options{ref} --pe";
if ( $interleaved == 1 ) {
    $mrfast_cmd = $mrfast_cmd . " --seq $options{fastq}";
}
else {
    my $fastq2 = defined $options{fastq2} ? $options{fastq2} : "$tmp_dir/$tmp_fn\_R2$split_number$fastq_suffix_name$gz_suf";
    if ( !-e $fastq2 ) { die "Unable to determine --fastq2 file. Please try using --fastq2. $fastq2\n"; }
    $mrfast_cmd = $mrfast_cmd . " --seq1 $options{fastq1} --seq2 $fastq2";
}
if ( $override_gzipped_status == 1 ) { $mrfast_cmd = $mrfast_cmd . " --seqcomp"; }
$mrfast_cmd = $mrfast_cmd . " --outcomp --discordant-vh --min $options{min} --max $options{max} -o $output";

my $append_lane_number = defined $options{append_lane_number} ? $options{append_lane_number} : "1";
my $auto_append        = defined $options{auto_append}        ? $options{auto_append}        : "1";
my $append_ids = ( defined $options{append_sample} or defined $options{append_rg} or defined $options{append_lib} ) ? "1" : $auto_append;

my $append_string;
if ( $append_ids == 1 ) {
    my $sample_id = defined $options{append_sample} ? $options{append_sample} : $split_number;
    if ( defined $options{append_sample} or $auto_append == 1 ) { $mrfast_cmd = $mrfast_cmd . " --sample $sample_id"; }

    my $rg;
    if   ( $append_lane_number == 1 ) { $rg = defined $options{append_rg} ? $options{append_rg} : $lane_number; }
    else                              { $rg = defined $options{append_rg} ? $options{append_rg} : $tmp_fn; }
    if ( defined $options{append_rg} or $auto_append == 1 ) { $mrfast_cmd = $mrfast_cmd . " --rg $rg"; }

    my $lib;
    if   ( $append_lane_number == 1 ) { $lib = defined $options{append_lib} ? $options{append_lib} : $tmp_fn; }
    else                              { $lib = defined $options{append_lib} ? $options{append_lib} : $options{ref}; }
    if ( defined $options{append_lib} or $auto_append == 1 ) { $mrfast_cmd = $mrfast_cmd . " --lib $lib"; }
}

run_cmd($mrfast_cmd);
if ( $clean_non_VH_files == 1 ) { run_cmd("rm $output.gz") }
print_complete( \%options );

sub help {
    die "\nThis script will run align fastq against a reference using mrfast.
	--fastq=	Interleaved Fastq
	--fastq1=	/path/to/R1.fastq.gz
	--fastq2=	Optional. If --fastq1 is used, the script will try to determine --fastq2. 
	--ref=
	--output_dir=
	--output_prefix=\n\n";
}
