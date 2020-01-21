## NEED TO ADD FASTQ DEDUPLICATION TO HERE TOO!

grep 'Total read' "${qc_dir}"/"${sample}".cutadapt.txt \
    | awk -v sample="${sample}" -v OFS='\t' '{print sample, "Total reads", $5}' \
    | tr -d ','
grep 'short' "${qc_dir}"/"${sample}".cutadapt.txt \
    | awk -v sample="${sample}" -v OFS='\t' '{print sample, "Reads < 20bp", $6}' \
    | tr -d ','

total_pairs=$(grep -m 1 'reads; of these:' "${qc_dir}"/"${sample}".alignment_stats.txt \
                | awk '{print $1}')
both_pair_unmapped=$(samtools view -c -f 12 -@ "${threads}" \
                        "${data_dir}"/"${sample}".fixmate.bam)
printf '%s\tBoth unmapped\t%s\n' "${sample}" "${both_pair_unmapped}"

r1_map_r2_unmap=$(samtools view -c -f 72 -F 4 -@ "${threads}" \
                    "${data_dir}"/"${sample}".fixmate.bam)
printf '%s\tOnly R1 mapped\t%s\n' "${sample}" "${r1_map_r2_unmap}"

r2_map_r1_unmap=$(samtools view -c -f 136 -F 4 -@ "${threads}" \
                    "${data_dir}"/"${sample}".fixmate.bam)
printf '%s\tOnly R2 mapped\t%s\n' "${sample}" "${r2_map_r1_unmap}"

unmapped_pairs=$(( both_pair_unmapped + r1_map_r2_ummap + r2_map_r1_unmap ))
duplicate_pairs=$(grep -m 1 'DUPLICATE' "${qc_dir}"/"${sample}".dedup_stats.txt \
                    | awk '{print $3}')
printf '%s\tDuplicates\t%s\n' "${sample}" "${duplicate_pairs}"

total_kept_pairs=$(samtools view -c -f 64 -@ "${threads}" \
                    "${data_dir}"/"${sample}".proc.bam)
qual_filtered_pairs=$(( total_pairs - total_kept_pairs - duplicate_pairs - unmapped_pairs ))
printf '%s\tMAPQ < 15\t%s\n' "${sample}" "${qual_filtered_pairs}"

