################################################
#File Name: circ_ratio.pl
#Author: Wanqiu Ding 
#Mail: wanqiuding@163.com
#Created Time: Wed Nov 25 14:48:56 2015
################################################

#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use vars qw($circ $junc $out);

my ($circ,$junc,$out);
GetOptions(
				"circ|c:s"		=>		\$circ,
				"junc|j:s"		=>		\$junc,
				"out|o:s"		=>		\$out,
);

open CIRC,"$circ" or die "can't open the file!";
open JUNC,"$junc" or die "can't open the file!";
open AMONG,">among.txt" or die "can't open the file!";
open OUT,">$out" or die "can't open the new file!";

my (%juncleft,%juncright);
while(<JUNC>){
	chomp;
	my ($chr,$chrStart,$chrEnd,$id,$numRead,$strand,$start,$end,$rbm,$num,$size,$block)=split('\t',$_);
	my @blockSize = split(',',$size);
	my @blockStart = split(',',$block);
	my $end5 = $blockSize[0] + $chrStart;
	my $end3 = $chrStart + $blockStart[-1];
	$juncleft{"$chr:$end5"}->{"$chr:$end3"}=$numRead;
	$juncright{"$chr:$end3"}->{"$chr:$end5"}=$numRead;
	print AMONG "$chr\t$end5\t$end3\t$id\t$numRead\t$strand\n";
	}
=cut
foreach my $keys1 (keys %junc){
	my $hash2 = $junc{$key1};
	foreach my $key2 (sort{$hash2->{$b}<=>$hash2->{$a}}keys %$hash2){
		print $key1."\t".$key2."\t".$hash2->{$key2}."\n";
		}
		}
=cut
while(<CIRC>){
	chomp;
	my @line = split "\t",$_;
	my ($right,$left)=(0,0);
	if(exists $juncleft{"$line[0]:$line[2]"}){
		my $hash2 = $juncleft{"$line[0]:$line[2]"};
		foreach my $key2 (sort{$hash2->{$b}<=>$hash2->{$a}}keys %$hash2)
		{
			$right+=$hash2->{$key2};
		}
	}else{
		$right=0;
	}
	if(exists $juncright{"$line[0]:$line[1]"}){
		my $hash2 = $juncright{"$line[0]:$line[1]"};
		foreach my $key2 (sort{$hash2->{$b}<=>$hash2->{$a}}keys %$hash2){
			$left+=$hash2->{$key2};
		}
	}else{
		$left=0;
	}
	print OUT join "\t",(@line,$left,$right);
	print OUT "\n";
	}
