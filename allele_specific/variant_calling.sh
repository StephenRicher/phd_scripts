#!/usr/bin/env bash

shopt -s extglob

print_usage() {
  echo "Usage: ..."
}

# Define default parameters that may be optionally specified
out=".";
threads=1

while getopts 's:g:r:h:v:o:t:' flag; do
  case "${flag}" in
    s) sample="${OPTARG}" ;;
    g) genome="${OPTARG}" ;;
    r) capture_regions="${OPTARG}" ;;
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

# If genome is not indexed then raise error.
if [ ! -f "${genome}".fai ] ; then
  echo "Genome must be indexed."; exit 
fi

# Merge coordinate sorted replicates. This is run first to ensure all threads are available from bcftools.
samtools merge -@ ${threads} \
  ${out}/${sample}.replicate_merge.bam "${@}"
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
    --output-file ${out}/${sample}.sorted.vcf.gz \
    --output-type z

bcftools index ${out}/${sample}.sorted.vcf.gz
	
# Subset VCF by capture region and perform haplotype phasing.
N=2
while IFS=$'\t' read -r chr start end region; do
  ((i=i%N)); ((i++==0)) && wait
  {
  bcftools view \
    --output-type v \
    --regions ${chr}:$((${start}+1))-${end} \
    ${out}/${sample}.sorted.vcf.gz \
    > ${out}/${sample}_${region}.vcf
  java -Xmx30g -jar ${hapcompass} \
    --bam ${out}/${sample}.replicate_merge.bam \
    --vcf ${out}/${sample}_${region}.vcf \
    --output ${out}/${sample}_${region}_hapcompass
  } &
done <"${capture_regions}"

MWER_solution_all="${out}/${sample}_all_hapcompass_MWER_solution.txt"
rm ${MWER_solution_all}
while IFS=  read -r -d $'\0' file; do
  # Get line number of header of top scoring section
  top_score_line=$(grep -n BLOCK ${file} | sort -k 6 -nr | cut -f 1 -d ':' | head -n 1)
  # Print out only the top scoring block
  awk -v n=${top_score_line} 'NR<n {next} NR==n {print;next} /^BLOCK/ {exit} {print}' ${file} >> ${MWER_solution_all}
done < <(find ${out} -type f -name "${sample}_*_hapcompass_MWER_solution.txt" -print0)

# Convert solution to VCF for genome masking
java -jar ${hc2vcf} ${MWER_solution_all_regions} <(gunzip -c ${out}/${sample}.sorted.vcf.gz) 2 true
	
# Convert VCF to SNPsplit friendly output
awk -v OFS=$'\t' 'substr($10,1,3)=="0|1" {print $3, $1, $2, 1, $4"/"$5} substr($10,1,3)=="1|0" {print $3, $1, $2, 1, $5"/"$4}' \
  ${MWER_solution_all_regions}.vcf \
  > ${out}/${sample}_MWER_solution_snpsplit.txt

## Masked reference genome at SNPs and build bowtie2 index ##

#genome_rmpath="${genome##*/}"
#masked_genome="${out}"/"${genome_rmpath%.fa?(.gz)}"_${sample}_masked.fa.gz
#masked_genome_index="${out}"/GRCh38_"${sample}"

#bedtools maskfasta \
#  -fullHeader \
#  -fi <(zcat -f "${genome}") \
#  -bed ${sample}_all_regions_hapcompass_MWER_solution.txt.vcf \
#  -fo /dev/stdout \
#| gzip > "${masked_genome}"

#bowtie2-build \
#  --threads "${threads}" \
#  "${masked_genome}" \
#  "${mask_genome_index}" 
