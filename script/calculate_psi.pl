#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

my $file = '';
my $key = '';
my $junc = '';
my $outfile = '';

GetOptions ('file=s'=>\$file, 'key=s'=>\$key,
			'junc=s'=>\$junc, 'outfile=s'=>\$outfile);

open (FILE,"$file") or die "can't open the uniqjuncpos file of ES!";
open (KEY,"$key") or die "can't open the junction key file!";
open (JUNC,"$junc") or die "can't open the total junction file!";
open (OUT,">$outfile") or die "can't open the new file!";

my (%hash,%juncHash);
while(<FILE>){
	chomp;
	next if($_ =~/^#/);
	my @field = split "\t",$_;
	my $file_key = join ('-',@field[1,2,5,16]);   #the key to indicate a circRNA
	$hash{$file_key}->{$field[3]} = {
		score => $field[4],
		sup_read => $field[12],
		total_junc_read => $field[15],
	};
}
while(<JUNC>){
	chomp;
	my @juncField = split "\t",$_;
	my @junc = split ":",$juncField[3];
	my $juncKey = join ":",(@junc[0,1,2]);
	$juncHash{$juncKey} = $junc[3];
}
	
while(<KEY>){
	chomp;
	my @field = split "\t",$_;
	my $join_junc = join "-",@field[6,7];
	my $new_junc = join ":",($field[0],$join_junc);
	my $line_key = join ('-',@field[1,2,5],$new_junc);
	my $existJunc = join ":",($field[0],$field[5],$join_junc);
	#my $line_junc_key = join (':',@field[11,4]); replaced by $field[7]
	my $up_junc;
	if($field[5] eq '+'){
		$up_junc = join ('-',$field[6],$field[1]);
	}else{
			$up_junc = join ('-',$field[2],$field[7]);
		}
	my $up_junc_key = join (':',$field[0],$field[5],$up_junc);
	if(exists $hash{$line_key}){
		my ($score,$sup_read,$total_read,$count) = (0,0,0,0);
		my $hash2 = $hash{$line_key};
		foreach my $key2 (sort{$hash2->{$b}<=>$hash2->{$a}}keys %$hash2)
		{
			$score+=$hash2->{$key2}->{score}/1000;
			$sup_read+=$hash2->{$key2}->{sup_read};
			$total_read+=$hash2->{$key2}->{total_junc_read};
			$count++;
		}
		print OUT join "\t",(@field,$score/$count,$sup_read/$count,$total_read/$count,"essup");
		print OUT "\n";
	}else{
		if(exists $juncHash{$existJunc}){
			print OUT join "\t",(@field,1,$juncHash{$existJunc},$juncHash{$existJunc},"juncsup");
			print OUT "\n";
		}elsif(exists $juncHash{$up_junc_key}){
			print OUT join "\t",(@field,0,0,$juncHash{$up_junc_key},"upjunc");
			print OUT "\n";
			}else{
				print OUT join "\t",(@field,"NA",0,0,"noexp");
				print OUT "\n";
			}
		}
	}
