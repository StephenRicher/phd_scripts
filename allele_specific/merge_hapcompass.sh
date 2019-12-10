#!/usr/bin/env bash

main() {
    local samples=( 'HB2_CL4' 'HB2_WT' 'MCF7' )
    data_dir="/home/u/sr467/scratch/projects/hic-01/allele_specific"

    for sample in "${samples[@]}"; do
        cat "${data_dir}"/"${sample}"*snpsplit* \
            > "${data_dir}"/"${sample}"_all_MWER_solution_snpsplit.txt

        for file in "${sample}"*MWER*solution*vcf; do
            bcftools view -O z "${file}" > "${file}".gz
            bcftools index "${file}".gz
        done

        bcftools concat --allow-overlaps  \
                "${data_dir}"/"${sample}"*MWER*solution*vcf.gz \
            > "${data_dir}"/"${sample}"-MWER_solution_all.vcf

    done

}

main
