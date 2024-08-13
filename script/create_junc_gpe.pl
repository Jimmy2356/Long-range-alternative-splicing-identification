#!/usr/bin/perl -w


#The purpose of this script is to generate all possible exon-exon junction annotations for each annotated gene
#the gpe annotation of junctions can be used for mapping step

use strict;
use Getopt::Long;
use File::Basename;

#Initialize command line options
# my $verbose = '';
my $outdir = '';
#my $logfile = '';
my $annotation = '';
#my $junc_gpe = '';

GetOptions ('outdir=s'=>\$outdir,'annotation=s'=>\$annotation);

my $junc_gpe = "$outdir"."exonJunctions".".gpe";
my $junc_tophat = "$outdir"."exonJunctions_"."list_juncs";

#open the gpe file
open (ANNO, "$annotation") || die "\nCould not open the annotation file:$annotation";
open (JUNC_GPE,">$junc_gpe") || die "\nCould not open the new file:$junc_gpe";
open (JUNC_TOPHAT,">$junc_tophat") || die "\nCould not open the new file:$junc_tophat";

while(<ANNO>){
	chomp;
	my @field = split "\t",$_;
	if($field[7] ==0){
		next();
		}
	my @exonS = split ',',$field[8];
	my @exonE = split ',',$field[9];	
	foreach my $exon1_count (0..($field[7]-2)){
		foreach my $exon2_count(($exon1_count+1)..($field[7]-1)){
			my $exon_Start = join ",",($exonS[$exon1_count],$exonS[$exon2_count]);
			my $exon_End = join ",",($exonE[$exon1_count],$exonE[$exon2_count]);
			print JUNC_GPE join "\t",("JUNC",$field[1],$field[2],$exonS[$exon1_count],$exonE[$exon2_count],$exonS[$exon1_count],$exonE[$exon2_count],2,$exon_Start,$exon_End);
			print JUNC_GPE "\n";
			print JUNC_TOPHAT join "\t",($field[1],$exonE[$exon1_count]-1,$exonS[$exon2_count],$field[2]);
			print JUNC_TOPHAT "\n";
			}
			}
}
		
		
		
		
		
		
		
		