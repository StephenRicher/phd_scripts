#!/usr/bin/env bash

genome=/home/stephen/x_db/DBuck/s_richer/genomes/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
dir=/home/stephen/x_db/DBuck/s_richer/genomes/GRCh38/ctcf

for ctcf in "${dir}"/MCF-7_CTCF*.txt; do
    ./extract_ctcf_direction.sh "${genome}" "${ctcf}" \
    > "${ctcf/.txt/.bed}"
done

# Merge CTCF replicates
cat "${dir}"/MCF-7_CTCF_*.bed \
| bedtools sort \
| bedtools merge -s -c 4,5,6 -o count,median,distinct \
> "${dir}"/MCF-7_CTCF-merged.bed
