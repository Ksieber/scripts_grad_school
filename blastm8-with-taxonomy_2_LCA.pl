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
	'input|i=s',
	'output|o=s',
	'lineage_col|c=i',
	'help|h',
	);




if ($options{help}){die "\nHELP: This script will take a file that is filtered for best hits from a blast output with taxonomy info and calculate LCA for each read.\n
   --input= blast.m8 report, filtered for best hits with lgt_filter_best_hit.pl to add taxon info.  
   --output= Name the output.
   --lineage_col= Column with lineage info. 
   --family_LCA : Will only the LCA at family lvl higher.\n"};

my $output_checker;
if($options{output}){open (OUTPUT, ">", $options{output}) || die "Can't open output because: $!\n"; $output_checker=$options{output};};

my %calculated_lca;
open (INPUT, "<", $options{input}) || die "ERROR: Couldn't open: $options{input} because $!\n";
while (<INPUT>){
    chomp;
    my @fields = split(/\t/, $_);
    my $read_id = $fields[0];
    my $lineage = $options{lineage_col} ? $options{lineage_col} : $fields[14];
    if ($calculated_lca{$read_id} && ($last_read_id eq $read_id)){
		my $lca = $calculated_lca{$read_id} ? $calculated_lca{$read_id} : $lineage;
		$lca = &find_lca($lineage, $lca);
		$calculated_lca{$read_id} = $lca;
    }
    if ($last_read_id ne $read_id){
    	if ($output_checker){ 
    		print OUTPUT "$last_read_id\t$calculated_lca{$last_read_id}\n";
		} else {
			print STDOUT "$last_read_id\t$calculated_lca{$last_read_id}\n";
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

# Print the final LCA. 
if ($output_checker){ 
    print OUTPUT "$last_read_id\t$calculated_lca{$last_read_id}\n";
} else {
	print STDOUT "$last_read_id\t$calculated_lca{$last_read_id}\n";
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

