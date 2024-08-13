#! /bin/bash
set -o errexit
### jimj
### Tue Oct 24 21:07:41 CST 2023
### step1: calculate junc from sam2jun.pl,
### step2: calculate scan file from scans.py
### step3: uniq circRNA file,calculate psi


### input
refdir='/home/user/data3/jimj/circRNA/circRNA_revision/simulation/ref'
scriptdir='/home/user/data3/jimj/circRNA/circRNA_revision/scripts/bin'


#eg.../A549/mapping/SRR2962431_tophat_fusion/circ.bed6+
#eg.../A549/mapping/SRR2962431_tophat_out/uniq.sorted.bam
export PERL5LIB=/home/user/data3/jimj/bin/:$PERL5LIB

### para
bam=/home/user/data3/jimj/circRNA/circRNA_revision/simulation/result/mixdata/mixdata_tophat_out/uniq.sorted.bam
circ=/home/user/data3/jimj/circRNA/circRNA_revision/simulation/result/mixdata/mixdata_tophat_fusion/circ.bed6+
prefix=mixdata
outdir=/home/user/data3/jimj/circRNA/circRNA_revision/simulation/result/linear/$prefix


	mkdir -p $outdir/junc
	mkdir -p $outdir/scans
	mkdir -p $outdir/psi
	mkdir -p $outdir/circRatio
	juncDir=$outdir/junc
	psiDir=$outdir/psi
	circDir=$outdir/circRatio
	scanDir=$outdir/scans
	
	
	## get junc
	perl $scriptdir/sam2junc.jimj.pl $bam > $juncDir/${prefix}_junc.bed12 2>$juncDir/${prefix}_junc.log 

	## get scans
	python $scriptdir/scan_es.py -j $juncDir/${prefix}_junc.bed12 \
	-g /home/user/data3/jimj/database/gpe/hg19.refGene.gpe -b \
	> $scanDir/${prefix}.scan_es_uniq.out
	2> $scanDir/${prefix}.scan_es_uniq.log

	python $scriptdir/es_freq_add.py -i $scanDir/${prefix}.scan_es_uniq.out \
	 -o $scanDir/${prefix}.scan_es_uniq.bed12 -g $scanDir/${prefix}.scan_es_uniq.group

	## get uniq ratio of circRNA
	#circ=$indir/$sample/mapping/${prefix}_tophat_fusion/circ.bed6+
	junc=$juncDir/${prefix}_junc.bed12
	scans=$scanDir/${prefix}.scan_es_uniq.bed12
	#cd $circDir
	sh $scriptdir/addRatioPsi.jmj.sh -c $circ -j $junc -s $scans -n $circDir/$prefix
	
	
wait 
echo "all is done"


