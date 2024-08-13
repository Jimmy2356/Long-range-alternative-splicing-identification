#!/bin/bash

## input
refdir='/rd1/user/jimj/circ_revision_20240726/compare_method/realdata/ref'
scriptdir='/rd1/user/jimj/circ_revision_20240726/compare_method/scripts/bin'

# para
fqfile=$1
fq2=$2
prefix=$3
outdir=$4
thread=$5

### mapping 
# under conda env py27

### unstrand
mkdir -p $outdir

        tophat \
        --read-mismatches 8 \
        --read-gap-length 3 \
        --read-edit-dist 8 \
        --read-realign-edit-dist 0 \
        --output-dir $outdir/${prefix}_tophat_out \
        --min-anchor-length 8 \
        --num-threads $thread \
        --no-coverage-search --microexon-search \
        --library-type fr-unstranded \
        --segment-mismatches 2 --segment-length 25 \
        -j $refdir/uniq_exonJunctions_list_juncs \
        $refdir/bowtie2/hg38 \
        $fqfile $fq2 \
        2>$outdir/${prefix}_tophat_out.log && \

        bamToFastq -i $outdir/${prefix}_tophat_out/unmapped.bam -fq $outdir/${prefix}_tophat_out/unmapped.fastq && \

        samtools view -@ $thread -bu -q 50 $outdir/${prefix}_tophat_out/accepted_hits.bam | samtools sort - -@ $thread -m 1G -o $outdir/${prefix}_tophat_out/uniq.sorted.bam \
        2>$outdir/${prefix}_tophat_out/samtoolsSort.log &&\

        samtools index -@ $thread $outdir/${prefix}_tophat_out/uniq.sorted.bam &&\

        tophat -o $outdir/${prefix}_tophat_fusion -p $thread \
        --fusion-search \
        --keep-fasta-order --bowtie1 --no-coverage-search \
        --library-type fr-unstranded \
        $refdir/bowtie1/hg38 \
        $outdir/${prefix}_tophat_out/unmapped.fastq \
        2>$outdir/${prefix}_tophat_fusion.log &&\

        perl $scriptdir/annotationhash_fusionJunc.pl \
        -f $outdir/${prefix}_tophat_fusion/accepted_hits.bam \
        -g $refdir/hg38.fa \
        -r $refdir/Homo_sapiens.GRCh38.103.gpe \
        -p $outdir/${prefix}_tophat_out/uniq.sorted.bam \
        -o $outdir/${prefix}_tophat_fusion/circ.bed6+ \
	2> $outdir/${prefix}_tophat_fusion/annotationhash_fusionJunc.log



wait
echo "done"




