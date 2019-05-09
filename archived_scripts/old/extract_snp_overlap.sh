#!/usr/bin/env bash

process_hicup() {

## Split interleaved BAM, sort and index.

# Remove file path and all extensions including the R1_4 part of filename.
rmext=${1%%_R1_4*}; sample=${rmext##*/}

for read_pair in R1 R2; do
{
        if [ "${read_pair}" == R1 ]; then flag="0x40"; else flag="0x80"; fi
	samtools view -hb -f ${flag} ${1} > ${sample}.${read_pair}.hicup.bam
	samtools sort ${sample}.${read_pair}.hicup.bam > ${sample}.${read_pair}.hicup.sorted.bam
	samtools index ${sample}.${read_pair}.hicup.sorted.bam
} &
done
}

extract_region() {

region=${2}; rmext=${1%%_R1_4*}; sample=${rmext##*/}

# Check if process_hicup output files are present. If not run to split, sort and index.
if [ ! -f ${sample}.R1.hicup.bam ] || [ ! -f ${sample}.R2.hicup.bam ] ; then process_hicup ${1}; fi

for read_pair in R1 R2; do
{
	samtools view -L ${region}_snp.bed ${sample}.${read_pair}.hicup.sorted.bam | cut -f 1 > ${sample}.${read_pair}.${region}.snp_overlap.txt
} &
done; wait

# Find read IDs overlapping the specified region that are common to both files and write to file.
comm -12 <(sort ${sample}.R1.${region}.snp_overlap.txt) <(sort ${sample}.R2.${region}.snp_overlap.txt) > ${sample}.both.${region}.snp_overlap.txt
					
for read_pair in R1 R2; do
{	
	samtools view ${sample}.${read_pair}.hicup.bam | LC_ALL=C grep -w -F -f ${sample}.both.${region}.snp_overlap.txt |
	cat <(samtools view -H ${sample}.${read_pair}.hicup.bam) - |
	samtools view -Sb > ${sample}.${read_pair}.${region}.both_snp.bam
} &	
done; wait


# Remove intermediate files
rm ${sample}.R1.${region}.snp_overlap.txt ${sample}.R2.${region}.snp_overlap.txt \
   ${sample}.both.${region}.snp_overlap.txt
}

snp=$1; region=$2; chr=$3; start=$4; end=$5

awk -v OFS='\t' -v chr=${chr} -v start=${start} -v end=${end} '$1 == chr && $3 >= start && $3 <= end {print}' ${snp} > ${region}_snp.bed

# Shift	to ensure the remaining	arguments are all the hicup files
shift 5

for bam in "${@}"; do
  	extract_region ${bam} ${region} ${region}_snp.bed &
done; wait
