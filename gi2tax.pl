#!/usr/bin/perl

=head1 NAME

gi2tax.pl

=head1 SYNOPSIS

USAGE:  gi2tax.pl --input /path/to/gis.list [--help]

=head1 OPTIONS

=over 10

=item B<--input>

List file with gi's in the format below:
gi|123456|ref|NC_012345|

=item B<--help>

Display the pod2usage page for this utility

=back

=head1 DESCRIPTION

gi2tax.pl - Reports taxonomic information for a list of gi's.

=cut

use FindBin;
use lib $FindBin::Bin;
use strict;
use GiTaxon;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Basename;
$|++;

my %options = ();
my $results = GetOptions (\%options,
              'input=s',
              'tmp_dir:s',
              'ncbitax=s',
              'gitax=s',
              'dbhost=s',
              'taxondb=s',
              'taxoncoll=s',
              'idx_dir=s',
              'taxon_dir=s',
              'help|h') || pod2usage();

# display documentation
if( $options{'help'}){
    pod2usage( {-exitval=>0, -verbose => 1, -output => \*STDOUT} );
}
if(!$options{'input'}) {
    pod2usage({-exitval=>0, -verbose => 1, -output => \*STDOUT});
}

my $gi2taxon = &get_gi2taxon();

open IN, "<$options{input}" or die $!;
print STDERR "About to parse the file\n";
while (my $line = <IN>) {
    chomp($line);
	my @fields = split(/\|/, $line);
    
    my $gi = $fields[1]; 
    my $tax = $gi2taxon->getTaxon($gi);

    my @newfields = ($line);

    if($tax->{'taxon_id'}) {
        push(@newfields,$tax->{'taxon_id'});
        if($tax->{'name'}) {
            push(@newfields,$tax->{'name'},$tax->{'lineage'});
        }
    }
    else {
        push(@newfields,('',''));

    }
    print join("\t",@newfields);
    print "\n";
}

sub get_gi2taxon {
    my $ncbitax = $options{ncbitax} ? $options{ncbitax} : '/local/db/repository/ncbi/blast/20110831_122720/taxonomy/taxdump/';
    my $gi2tax = $options{gitax} ? $options{gitax} : '/local/db/repository/ncbi/blast/20110831_122720/taxonomy/gi_taxid_nucl.dmp';
    my $dbhost = $options{dbhost} ? $options{dbhost} : 'skimaro-lx.igs.umaryland.edu';
    my $taxondb = $options{taxondb} ? $options{taxondb} : 'gi2taxon';

    my $idx_dir = $options{idx_dir};

    if(!$idx_dir && -e "$ncbitax/names") {
        $idx_dir = $ncbitax;
    }
    elsif(!$idx_dir) {
        $idx_dir='/local/projects/HLGT/tax_idx/20110831/';
}

    my $gi2taxon = GiTaxon->new(
         {'taxon_dir' => $options{taxon_dir},
         'idx_dir' => $idx_dir,
         'host' => $dbhost,
        });

    return $gi2taxon;
}

