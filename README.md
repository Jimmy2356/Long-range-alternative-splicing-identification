# Long-range-alternative-splicing-identification
## Requirement
### Software
Recommond to install these software using conda, ensure this tools are in the environment path, else please replace the full path of the tools in the step1_mapping.sh
1. python 2.7
2. Tophat v2.1.0
3. bowtie v2.2.5
4. bowtie version v1.1.2
5. bedtools bamtofastq
6. perl bioperl (Bio/DB/Fasta.pm)
7. samtools v1.17
### File
1. Junction file for mapping with tophat. This step is to generate junction database for tophat -j option while mapping
```
perl script/create_junc_gpe.pl --outdir=ref/ --annotation=ref/Homo_sapiens.GRCh38.103.gpe 
less ref/exonJunctions_list_juncs | sort | uniq > ref/uniq_exonJunctions_list_juncs 
```
2. index of bowtie1, bowtie2, annotation of genome(.gpe)
## How to use
1. Mapping
step1_mapping.sh is used to mapping rRNA- RNA seq pair-end reads and output circRNA list (circ.bed6+) in *tophat_fusion directory. 
Before running the script, $refdir (include bowtie, bowtie2 index, genome fasta, .gpe annotation) and $scriptdir (include the script supply in this site) need to replaced with your path.
```
sh -x step1_mapping.sh $fq1 $fq2 $prefix $outdir > $log
# eg. sh -x realdata_mapping_nortrim.sh fastq/SRR17235470_1.fastq.gz fastq/SRR17235470_2.fastq.gz SRR17235470 results_nortrim/SRR17235470 > results_nortrim/SRR17235470.log 2>results_nortrim/SRR17235470.err &
```
2. identify long range splicing
step2_quantify_circRNA_and_LR.sh is used to identify and quantify circRNA and long range splicing transcript, output *.uniq.ratio.psi.bed6+ file in circRatio directory.
Before running the script, addRatioPsi.sh, the variable in "###input" and "###para" need to be replaced with your path.
```
sh -x step2_quantify_circRNA_and_LR.sh
```
The output file will look like below, circ_ratio and psi are the expression of the circRNA and corresponding long range splicing transcript
#chr    start   end     circ_id backsplice_supp strand  left_junc       right_junc      max_juncNum     circ_ratio      psi     flank_juncNum   total_juncNum   label
chr10   100013309       100016704       circ_1  11      -       100012225       100017406       70      0.135802        0.102   8       78      essup
chr10   100013309       100017921       circ_2  24      -       100012225       100018765       81      0.228571        0.19    19      100     essup
chr10   100013309       100020884       circ_3  24      -       100012225       100021791       39      0.380952        0.264   14      53      essup
chr10   100152188       100152840       circ_4  24      -       100150839       100154952       13      0.648649        0.566   17      30      essup
chr10   101439017       101439632       circ_5  7       +       101421385       101445548       46      0.132075        0.22    13      59      essup
chr10   101507013       101507147       circ_6  5       +       101503829       101510125       20      0.2     0.13    3       23      essup


