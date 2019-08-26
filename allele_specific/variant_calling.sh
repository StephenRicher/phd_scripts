#!/usr/bin/env bash

shopt -s extglob

print_usage() {
  echo "Usage: ..."
}

# Define default parameters that may be optionally specified
out=".";
threads=1

while getopts 's:g:r:f:b:h:v:o:t:' flag; do
  case "${flag}" in
    s) sample="${OPTARG}" ;;
    g) genome="${OPTARG}" ;;
    r) capture_regions="${OPTARG}" ;;
    f) fastqs+=("$OPTARG") ;;
    b) bams+=("$OPTARG") ;;
    h) hapcompass="${OPTARG}" ;;
    n) snpsplit="${OPTARG}" ;;
    v) hc2vcf="${OPTARG}" ;;
    o) out="${OPTARG}" ;;
    t) threads="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

# Shift to access all other unparsed arguments (all bam files to merge)
shift "$((OPTIND-1))"

if [ ! "${#fastqs[@]}" -eq "2" ]; then
  >&2 echo "Expecting two fastq files for -f."
  exit 1
fi

if [ ! "${#bams[@]}" -eq "2" ]; then
  >&2 echo "Expecting two bam files for -b"
  exit 1
fi

# If genome is not indexed then raise error.
if [ ! -f "${genome}".fai ] ; then
  echo "Genome must be indexed."; exit 
fi

# Create directory to store VCF output.
vcfs="${out}/vcfs"
mkdir ${vcfs}

# Merge coordinate sorted replicates. This is run first to ensure all threads are available from bcftools.
samtools merge -@ ${threads} \
  ${out}/${sample}.replicate_merge.bam "${bams[@]}"
samtools index -@ ${threads} \
  ${out}/${sample}.replicate_merge.bam

# Call variants from merged sample file and quality filter.
bcftools mpileup \
  -q 15 --ignore-RG --count-orphans \
  --max-depth 100000 \
  --output-type u
  -f "${genome}" \
  --regions-file ${capture_regions} \
  "${out}"/"${sample}".replicate_merge.bam \
| bcftools call \
    --skip-variants indels \
    --multiallelic-caller \
    --variants-only \
    --output-type u \
| bcftools view \
    -i '%QUAL>=20' \
    --output-type u \
| bcftools sort \
    --output-file ${vcfs}/${sample}.sorted.vcf.gz \
    --output-type z

bcftools index ${vcfs}/${sample}.sorted.vcf.gz
	
# Subset VCF by capture region and perform haplotype phasing.
N=2
while IFS=$'\t' read -r chr start end region; do
  ((i=i%N)); ((i++==0)) && wait
  {
  bcftools view \
    --output-type v \
    --regions ${chr}:$((${start}+1))-${end} \
    ${vcfs}/${sample}.sorted.vcf.gz \
    > ${vcfs}/${sample}_${region}.vcf
  java -Xmx30g -jar ${hapcompass} \
    --bam ${out}/${sample}.replicate_merge.bam \
    --vcf ${vcfs}/${sample}_${region}.vcf \
    --output ${vcfs}/${sample}_${region}_hapcompass
  } &
done <"${capture_regions}"

MWER_solution_all="${vcfs}/${sample}_all_hapcompass_MWER_solution.txt"
rm ${MWER_solution_all}
while IFS=  read -r -d $'\0' file; do
  # Get line number of header of top scoring section
  top_score_line=$(grep -n BLOCK ${file} | sort -k 6 -nr | cut -f 1 -d ':' | head -n 1)
  # Print out only the top scoring block
  awk -v n=${top_score_line} 'NR<n {next} NR==n {print;next} /^BLOCK/ {exit} {print}' ${file} >> ${MWER_solution_all}
done < <(find ${vcfs} -type f -name "${sample}_*_hapcompass_MWER_solution.txt" -print0)

# Convert solution to VCF for genome masking
java -jar ${hc2vcf} ${MWER_solution_all_regions} <(gunzip -c ${vcfs}/${sample}.sorted.vcf.gz) 2 true
	
# Convert VCF to SNPsplit friendly output
awk -v OFS=$'\t' 'substr($10,1,3)=="0|1" {print $3, $1, $2, 1, $4"/"$5} substr($10,1,3)=="1|0" {print $3, $1, $2, 1, $5"/"$4}' \
  ${MWER_solution_all_regions}.vcf \
  > ${vcfs}/${sample}_MWER_solution_snpsplit.txt

## Masked reference genome at SNPs and build bowtie2 index ##

genome_rmpath="${genome##*/}"
masked_genome_dir="${out}"/masked_genomes
masked_genome="${masked_genome}"/"${genome_rmpath%.fa?(.gz)}"_${sample}_masked.fa.gz
masked_genome_index_dir="${out}"/masked_genomes/index
masked_genome_index="${masked_genome_index_dir}"/GRCh38_"${sample}"
mkdir -p "${mask_genome_dir}" "${mask_genome_index_dir}"

bedtools maskfasta \
  -fullHeader \
  -fi <(zcat -f "${genome}") \
  -bed ${sample}_all_regions_hapcompass_MWER_solution.txt.vcf \
  -fo /dev/stdout \
| gzip > "${masked_genome}"

bowtie2-build \
  --threads "${threads}" \
  "${masked_genome}" \
  "${mask_genome_index}" 

## Map and filter ##

qc="${out}"/qc
data_dir="${out}"/data
mkdir -p "${qc}" "${data_dir}"

alignment_stats="${qc}"/alignment_stats.txt

if [ ! -f "${alignment_stats}" ]; then
  printf 'sample\ttotal_mapped_pairs\tboth_pairs_unmapped\t' >> "${alignment_stats}"
  printf 'r1_map_r2_unmapped\tr2_map_r1_unmapped\tduplicate_pairs\t' >> "${alignment_stats}"
  printf 'quality_filtered_pairs\ttotal_pairs_retained\n' >> "${alignment_stats}"
fi

# File for intermediate logfile
intermediate="${data_dir}"/"${sample}".fixmate.bam

hictools map \
  --index "${mask_genome_index}"  \
  --sample "${sample}"
  --log "${qc}"/"${sample}".bowtie2.logfile \
  --intermediate "${intermediate}" \
  "${fastqs[@]}" \
  2> "${qc}"/"${sample}"_alignment_stats.txt \
| hictools deduplicate \
  --log "${qc}"/"${sample}".dedup.logfile \
  2> "${qc}"/"${sample}".dedup_stats.txt \
| hictools process \
  --digest "${genome_digest}" \
  --log "${qc}"/"${sample}".process.logfile \
  > "${data_dir}"/"${sample}".proc.bam

total_pairs=$(grep -m 1 'reads; of these:' "${sample}".bowtie2_stats.txt | awk '{print $1}')
both_pair_unmapped=$(samtools view -cf 12 "${intermediate}")
r1_map_r2_unmap=$(samtools view -c -f 72 -F 4 "${intermediate}")
r2_map_r1_unmap=$(samtools view -c -f 136 -F 4 "${intermediate}")
unmapped_pairs=$(( both_pair_unmapped + r1_map_r2_ummap + r2_map_r1_unmap ))
duplicate_pairs=$(grep -m 1 'DUPLICATE' "${sample}".dedup_stats.txt | awk '{print $3/2}')
total_kept_pairs=$(samtools view -cf 64 "${sample}".proc.bam)
qual_filtered_pairs=$(( total_pairs - total_kept_pairs - duplicate_pairs - unmapped_pairs ))

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${sample}" "${total_pairs}" "${both_pair_unmapped}" \
  "${r1_map_r2_unmap}" "${r2_map_r1_unmap}" \
  "${duplicate_pairs}" "${qual_filtered_pairs}" \
  "${total_kept_pairs}" >> "${alignment_stats}"

"${snpsplit}" --hic \
  --snp_file "${vcfs}"/${sample}_MWER_solution_snpsplit.txt \
  -o "${data_dir}" \
  "${data_dir}"/"${sample}".proc.bam

# Extract replicate number from sample name
replicate_num=$(tr '-' '\n' <<<"${sample}" | head -n 2 | tail -n 1)
samtools merge -n -@ "${threads}" \
  "${sample}"_G1-${replicate_num}.bam \
  "${data_dir}"/"${sample}".proc.G1_!(G2).bam
samtools merge -n -@ "${threads}" \
  "${sample}"_G2-${replicate_num}.bam \
  "${data_dir}"/"${sample}".proc.G2_!(G1).bam
