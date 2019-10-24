#!/usr/bin/env bash

sample="${1}"

allele_dir=/home/u/sr467/scratch/projects/HiC/allele_specific/hic_dir
qc="${allele_dir}"
data_dir="${allele_dir}"
unmasked_genome_digest="${allele_dir}"/GRCh38_Mbo1-digest.txt.gz
masked_genome_index="${allele_dir}"/GRCh38_HB2_WT
snpsplit=/home/u/sr467/scratch/projects/HiC/allele_specific/HB2_WT_MWER_solution_snpsplit.txt
hictools="/home/u/sr467/scratch/scripts/hictools/hictools/hictools.py"

alignment_stats="${qc}"/alignment_stats.txt
printf 'sample\ttotal_mapped_pairs\tboth_pairs_unmapped\t' >> "${alignment_stats}"
printf 'r1_map_r2_unmapped\tr2_map_r1_unmapped\tduplicate_pairs\t' >> "${alignment_stats}"
printf 'quality_filtered_pairs\ttotal_pairs_retained\n' >> "${alignment_stats}"

# File for intermediate logfile
intermediate="${allele_dir}"/"${sample}".fixmate.bam

"${hictools}" map \
  --index "${masked_genome_index}"  \
  --sample "${sample}" \
  --log "${qc}"/"${sample}".bowtie2.logfile \
  --intermediate "${intermediate}" \
  --sensitivity sensitive \
  -@ 16 \
  "${data_dir}"/"${sample}"-R[14]-trim-trunc.fq.gz \
  2> "${qc}"/"${sample}".alignment_stats.txt \
| "${hictools}" deduplicate \
  --log "${qc}"/"${sample}".dedup.logfile \
  -@ 16 \
  2> "${qc}"/"${sample}".dedup_stats.txt \
| "${hictools}" process \
  --digest "${unmasked_genome_digest}" \
  --log "${qc}"/"${sample}".process.logfile \
  > "${allele_dir}"/"${sample}".proc.bam

"${hictools}" filter \
  --min_ditag 100 --max_ditag 1000 \
  --min_inward 1000 \
  --log "${qc}"/"${sample}".filter.logfile \
  --sample "${sample}" \
  "${allele_dir}"/"${sample}".proc.bam \
  > "${allele_dir}"/"${sample}".filt.bam \
  2>> "${qc}"/filter_statistics.tsv

total_pairs=$(grep -m 1 'reads; of these:' "${qc}"/"${sample}".alignment_stats.txt | awk '{print $1}')
both_pair_unmapped=$(samtools view -cf 12 "${intermediate}")
r1_map_r2_unmap=$(samtools view -c -f 72 -F 4 "${intermediate}")
r2_map_r1_unmap=$(samtools view -c -f 136 -F 4 "${intermediate}")
unmapped_pairs=$(( both_pair_unmapped + r1_map_r2_ummap + r2_map_r1_unmap ))
duplicate_pairs=$(grep -m 1 'DUPLICATE' "${qc}"/"${sample}".dedup_stats.txt | awk '{print $3}')
total_kept_pairs=$(samtools view -cf 64 "${allele_dir}"/"${sample}".proc.bam)
qual_filtered_pairs=$(( total_pairs - total_kept_pairs - duplicate_pairs - unmapped_pairs ))

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${sample}" "${total_pairs}" "${both_pair_unmapped}" \
  "${r1_map_r2_unmap}" "${r2_map_r1_unmap}" \
  "${duplicate_pairs}" "${qual_filtered_pairs}" \
  "${total_kept_pairs}" >> "${alignment_stats}"

/home/u/sr467/scratch/scripts/phd_scripts/allele_specific/SNPsplit.pl \
  --snp_file "${snpsplit}" \
  --hic \
  "${allele_dir}"/"${sample}".filt.bam

for a in 1 2; do
  samtools merge -n "${sample/-*/_G${a}-${sample/*-/}}".bam \
    "${sample}"*G"${a}"_[GU]["${a}"A].bam
done
