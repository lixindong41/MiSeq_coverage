#!/usr/bin/env bash
# Pipeline for analysis of HIV MiSeq data
# bash coverage.sh -1 SRR961514_1.fastq -2 SRR961514_2.fastq -x HXB2.fasta -q 30

#--- Read command line args
while getopts 1:2:x:q: option
do
 case "${option}"
 in
 1) miseq1=${OPTARG};;
 2) miseq2=${OPTARG};;
 x) reference=${OPTARG};;
 q) cutoff=${OPTARG};;
 esac
done
cutoff=${cutoff:-0}
if [ $cutoff -gt 41 ]
then
	echo "q is reset to its maximum threshold: 41"
	cutoff=41	
fi

#--- Initial variables
base=$(dirname "$0")
data=${base}/rawdata
ref=${base}/ref/${reference/.fasta/}
result=${base}/result
miseq1=${miseq1/.fastq/}
miseq2=${miseq2/.fastq/}

echo "STARTING TIME: [`date`]"

#--- Trimming paired-end reads
#cutadapt -a ADAPTER_FWD \
#	-A ADAPTER_REV \
#	-u 15 \
#	-o ${result}/${miseq1}_trimed.fastq \
#	-p ${result}/${miseq2}_trimed.fastq \
#	-q ${cutoff} \
#${data}/${miseq1}.fastq $data/${miseq2}.fastq

cutadapt -q ${cutoff} -a ADAPTER_FWD \
	-u 15 \
	-o ${result}/${miseq1}_trimed.fastq ${data}/${miseq1}.fastq
cutadapt -q ${cutoff} -a ADAPTER_REV \
	-u 15 \
	-o ${result}/${miseq2}_trimed.fastq ${data}/${miseq2}.fastq

#--- Align reads to the genome
bowtie2 -x ${ref} \
	-1 $result/${miseq1}_trimed.fastq \
	-2 $result/${miseq2}_trimed.fastq \
	--very-sensitive \
	-S $result/${miseq1}.sam

#--- Convert file .sam to .bam
samtools view \
	-bS $result/${miseq1}.sam \
> $result/${miseq1}.bam

# Remove PCR duplicates
samtools rmdup -sS $result/${miseq1}.bam $result/${miseq1}_dedup.bam 

samtools sort $result/${miseq1}_dedup.bam \
	-o $result/${miseq1}_sorted.bam

samtools index $result/${miseq1}_sorted.bam


#--- Count coverage across the genome
#samtools depth -q ${cutoff} $result/${miseq1}_sorted.bam \
#| awk -F '\t' -v OFS='\t' '$3 > 5 { print $2, $3 } ' \
#> ${result}/q${cutoff}_coverage.tsv

# generate a fasta index
samtools faidx ${ref}.fasta 

# extract sequence names (chromosomes) and their lengths to a new file
cut -f 1,2 ${ref}.fasta.fai > ${ref}.sizes

bedtools makewindows \
	-g ${ref}.sizes \
	-w 1 \
> ${ref}.bed
bedtools coverage -counts \
	-b $result/${miseq1}_sorted.bam  \
	-a ${ref}.bed \
| awk -F '\t' -v OFS='\t' '$4 > 5 { print $3, $4 } ' \
> ${result}/q${cutoff}_coverage.tsv

#--- Get the nucleotide content across the genome
# create a BED file with windows of wished width
width=100
bedtools makewindows \
	-g ${ref}.sizes \
	-w ${width} \
> ${ref}_${width}bps.bed

bedtools nuc \
	-fi ${ref}.fasta \
	-bed ${ref}_${width}bps.bed \
| awk -F '\t' -v OFS='\t' 'FNR > 1 { print $1, $2, $3, $6, $7, $8, $9 } ' \
> ${ref}_nuc.txt

bedtools coverage -counts \
	-b $result/${miseq1}_sorted.bam \
	-a ${ref}_${width}bps.bed \
> $result/${miseq1}_counts.txt

#--- Combine nuc.txt and counts.txt
awk -F '\t' -v OFS='\t' \
'NR==1 { print "chr", "start", "end", "num_A", "num_C", "num_G", "num_T", "coverage"; next } \
NR==FNR{a[$1,$2,$3]=$0;next} ($1 SUBSEP $2 SUBSEP $3) in a { print a[$1,$2,$3], $4 }' \
${ref}_nuc.txt $result/${miseq1}_counts.txt \
>$result/${miseq1}_q${cutoff}_corr.txt

#--- Call R to plot coverage and correlation
echo "${result}/q${cutoff}_coverage.tsv​​"
Rscript ./resultplot.R ${result}/q${cutoff}_coverage.tsv​​ ${result}/${miseq1}_q${cutoff}_corr.txt

echo "END TIME: [`date`]"


