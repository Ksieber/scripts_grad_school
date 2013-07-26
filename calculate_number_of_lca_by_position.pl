#!/usr/local/bin/perl -w
use strict;

my $line_count = 1;
my %position_tracker;
my %lca_unique;


while(<>){
    chomp; 
    my ($position, $lcas_seen) = split(/\t/,$_);
    $lcas_seen=~s/,\s/_/g;
    $lcas_seen=~s/ /_/g;
    my @diff_lca=split(',',$lcas_seen);
    foreach my $uniq_lca(@diff_lca){
        my ($number_of_lca,$seen_lca) = split(/:/,$uniq_lca);
        my @lca_breakdown=split(/;/,$seen_lca);
        my $specific_lca = $lca_breakdown[-1] ? $lca_breakdown[-1] : "None";
        $position_tracker{$position}->{$specific_lca}=$number_of_lca;
    }
    $line_count++;
}

for my $position_keys (sort keys %position_tracker){
    for my $lca_keys (sort keys %{%position_tracker->{$position_keys}}){
        $lca_unique{$lca_keys}++;
    }
}

## Print Header to the output w/ Position \t each LCA found
print STDOUT "Position";
for my $each_lca (sort keys %lca_unique){
    print STDOUT "\t$each_lca"
}
print STDOUT "\n";

## Print actual output. Go through each position, print all the # LCAs found
for (my $print_position = 1; $print_position < $line_count; $print_position++){
    print STDOUT "$print_position";
    for my $each_lca (sort keys %lca_unique){
        if ($position_tracker{$print_position}->{$each_lca}){
            print STDOUT "\t$position_tracker{$print_position}->{$each_lca}";
        } else {
            print STDOUT "\t0";
        }
    }
    print STDOUT "\n";
}




















__END__
while(<>){
   chomp;
   my ($position, $lcas)=split(/\t/,$_);
   print STDOUT "$position";
   my @diff_lca=split(/,/,$lcas);
   foreach my $uniq_lca(@diff_lca){
      my ($number_of_lca,$lca)=split(/:/,$uniq_lca);
      my @split_lca=split(/;/,$lca);
      print STDOUT "\t$number_of_lca\t$split_lca[-1]";
   }
   print STDOUT "\n";
}
