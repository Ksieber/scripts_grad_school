#!/usr/local/bin/perl 
##This script will take blast results that have been parsed for best hit and have taxonomy info (field 14) and create a file for LCA of each read (must be unique).
##INPUT = _overall.out from /local/projects/ergatis/package-driley/bin/filter_lgt_best_hit
##OUTPUT is STDOUT unless specified w/ --output

use strict;
use File::Basename;
use List::Util qw[min max];
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
$| =1;

my %options;
my $results = GetOptions (\%options,
	'input=s',
	'output=s',
	'lgt',
	'family_LCA',
	'help|h',
	);

my $output_checker;
my $last_read_id = "space_holder";
my %calculated_lca;

if ($options{help}){die "\nHELP: This script will take a file that is filtered for best hits from a blast output with taxonomy info and calculate LCA for each read.\n
   --input= blast.m8 report, filtered for best hits with lgt_filter_best_hit.pl to add taxon info.  
   --output= Name the output.
   --family_LCA : Will only the LCA at family lvl higher.\n"};

if ($options{output}){open (OUTPUT, ">", $options{output}) || die "Can't open output because: $!\n"; $output_checker=$options{output};};

open (INPUT, "<", $options{input}) || die "ERROR: Couldn't open: $options{input} because $!\n";
while (<INPUT>){
    chomp;
    my @fields = split(/\t/, $_);
    my $read_id = $fields[0];
    my $lineage = $fields[14];
    if ($calculated_lca{$read_id} && ($last_read_id eq $read_id)){
	my $lca = $calculated_lca{$read_id} ? $calculated_lca{$read_id} : $lineage;
	$lca = &find_lca($lineage, $lca);
	$calculated_lca{$read_id} = $lca;
    }
    if ($last_read_id ne $read_id){
	if ($output_checker){ 
	    if ($last_read_id ne "space_holder"){
		if (!$options{family_LCA}){
		    print OUTPUT "$last_read_id\t$calculated_lca{$last_read_id}\n";
		}
		if ($options{family_LCA}){
		    if ($calculated_lca{$last_read_id}=~/;(\w*aceae);/){
			print OUTPUT "$last_read_id\t$1\n";
		    } else {
			my @lca_lineage = split(';',$calculated_lca{$last_read_id});
			my $z = @lca_lineage;
			print OUTPUT "$last_read_id\t$lca_lineage[$z-1]\n";
		    }
		}
	    }
	}
	elsif ($last_read_id ne "space_holder") {
	    if (!$options{family_LCA}){
		print OUTPUT "$last_read_id\t$calculated_lca{$last_read_id}\n";
	    } 
	    if ($options{family_LCA}){
		if ($calculated_lca{$last_read_id}=~/;(\w*aceae);/){
		    print OUTPUT "$last_read_id\t$1\n";
		} else {
		    my @lca_lineage = split(';',$calculated_lca{$last_read_id});
		    my $z = @lca_lineage;
		    print OUTPUT "$last_read_id\t$lca_lineage[$z-1]\n";
		}
	    } 
	}
	for (keys %calculated_lca){
	    delete $calculated_lca{$_};
	}
	my $lca = $calculated_lca{$read_id} ? $calculated_lca{$read_id} : $lineage;
	$lca = &find_lca($lineage, $lca);
	$calculated_lca{$read_id} = $lca;
	$last_read_id = $read_id;
    } 
}


sub find_lca {
  my @lineages = shift;
  my @lca = split(';', $_[0]);
  
  foreach my $l (@lineages) {
    my $newlca = [];
    my @lineage = split(';',$l);
    for( my $i = 0; $i < @lineage; $i++) {
      if($lca[$i] eq $lineage[$i]) {
	push(@$newlca, $lineage[$i]);
      }
      else {
	last;
      }   
    }
    @lca = @$newlca;
  }
  return join(';',@lca);
}


__END__

