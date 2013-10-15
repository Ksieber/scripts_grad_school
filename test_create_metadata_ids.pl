#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;

my $count = 10000;
my %patient_id;
my %bam_id;
print STDOUT join ("\t",('kbs_bam_id','tcga_bam_id','kbs_patient_id','tcga_patient_id','Tissue_Type','unkown_id_number'));
print STDOUT "\n";
 

while(<>){
	chomp;
	my ($patient,$type,$unknown,$bam)=(split /\//, $_)[4,5,6,7];
	if( !$patient_id{$patient} ){ $patient_id{$patient} = $count; $count+=1000;}
	$bam_id{$bam}=$patient_id{$patient};
	if($type=~/tumor/){$bam_id{$bam}+=100;}	
	if($type=~/rna_seq/){$bam_id{$bam}+=200;}
	if($type=~/normal/){$bam_id{$bam}+=300;}
	`mkdir -p CRC$patient_id{$patient}`;
	`ln -s $_ CRC$patient_id{$patient}/CRC$bam_id{$bam}\.bam`;
	print STDOUT join("\t",("CRC$bam_id{$bam}\.bam","$bam","$patient_id{$patient}","$patient","$type","$unknown"));
	print STDOUT "\n";
}