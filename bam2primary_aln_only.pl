#!/usr/bin/perl
use lib ( '/home/ksieber/scripts/', '/home/ksieber/perl5/lib/perl5/' );
use warnings;
use strict;
use File::Basename;
use print_call;
use mk_dir;
use run_cmd;
use setup_input;

if ( not @ARGV ) { &help; }
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions( \%options, 'input|i=s', 'output_prefix|p=s', 'output_dir|o=s', 'output|O=s', 'Qsub|q=i', 'help|?', 'sub_mail=s', )
    or die "Error: Unrecognized command line option. Please try again.\n";

if ( $options{help} ) { &help; }
$options{input} = defined $options{input} ? $options{input} : $ARGV[0];
if ( !$options{input} && !$options{input_list} ) { die "Error: Must pass an input file with --input, --input_list, or $ARGV[0]. Please try again.\n"; }
my $input = $options{input};

Qsub_script( \%options ) if ( defined $options{Qsub} and $options{Qsub} == 1 );

my ( $in_fn, $in_dir, $in_suf ) = fileparse( $input, "\.bam" );
my $out_name = defined $options{output_prefix} ? $options{output_prefix} : "$in_fn\_prime";
my $out_dir  = defined $options{output_dir}    ? $options{output_dir}    : $in_dir;
$out_dir =~ s/\/+$//;
mk_dir($out_dir);
my $stdout = ( defined $options{output} and $options{output} =~ /STDOUT/i ) ? "1" : 0;
my $out = ( defined $options{output} and $stdout != 1 ) ? $options{output} : "$out_dir\/$out_name\.bam";

if ( $stdout == 1 ) {
    open( my $out_fh, "|-", "samtools view -S - -u" )         or die "Error: Unable to open the output filehandle.\n";
    open( my $in_fh,  "-|", "samtools view -hF 2304 $input" ) or die "Error: Unable to open the input filehand.\n";
    while (<$in_fh>) { print $out_fh "$_"; }
    close $in_fh;
    close $out_fh;
}
else {
    run_cmd("samtools view -F 2304 $input -bo $out");
}

print_complete( \%options );

sub help {
    die
        "\nThis script will parse a bam file, removing all secondary alignments. It is quicker to run samtools directly (samtools view -F0x100 \$input -bo \$out), but this allows for easy output naming and Qsub. 
	--input|i=			Input bam file to parse. May also take input from ARGV.
	--output|O=			< /full/path/and/name/for/output.bam > 	[input_dir/input_name_prime.bam] 
		--output can be set to \'STDOUT\' for piping to other commands, outputs uncompressed BAM
	--output_dir|o=		< /path/for/output > 				[input_dir]
	--output_prefix=	< prefix.bam > 						[(input)_prime.bam]
	--Qsub=				< 0|1 > [0] 1= Submit the job the SGE. 
	--help|?\n\n";
}

__END__
my ( $out_header, $in_fh ) = open_bam($input);
my $out_fh;
if($stdout==1){
	open( $out_fh, "|-", "samtools view -S - -u") or die "Error: Unable to open the output filehandle.\n";
	print $out_fh "$out_header\n"; 
}
else {
	$out_fh = write_bam( $out, $out_header )
}

while ( my $read = read_bam($in_fh) ) {
    next if ( $read->{flag}->{secondary} );
    print $out_fh "$read->{line}\n";
}
