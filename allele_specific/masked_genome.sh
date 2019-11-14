#!/usr/bin/env bash

shopt -s extglob

print_usage() {
  echo "Usage: ..."
}

# Define default parameters that may be optionally specified
out=".";
threads=1

while getopts 's:g:v:o:t:' flag; do
  case "${flag}" in
    s) sample="${OPTARG}" ;;
    g) genome="${OPTARG}" ;;
    v) vcf="${OPTARG}" ;;
    o) out="${OPTARG%/}" ;;
    t) threads="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

# If genome is not indexed then index.
if [ ! -f "${genome}".fai ] ; then
  samtools faidx "${genome}"
fi

## Masked reference genome at SNPs and build bowtie2 index ##

genome_rmpath="${genome##*/}"
masked_genome="${out}"/"${genome_rmpath%.fa?(.gz)}"_"${sample}"_masked.fa

bedtools maskfasta \
  -fullHeader \
  -fi <(zcat -f "${genome}") \
  -bed "${vcf}" \
  -fo "${masked_genome}"

bowtie2-build \
  --threads "${threads}" \
  "${masked_genome}" \
  "${out}"/GRCh38_"${sample}" 
