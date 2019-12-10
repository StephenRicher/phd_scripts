#!/usr/bin/env bash

main() {
    local samples=( 'HB2_CL4' 'HB2_WT' 'MCF7' )
    data_dir="/home/u/sr467/scratch/projects/hic-01/allele_specific"

    for sample in "${samples[@]}"; do
        cat "${data_dir}"/"${sample}"*snpsplit* \
            > "${data_dir}"/"${sample}"_all_MWER_solution_snpsplit.txt

        for file in "${sample}"*MWER*solution*vcf; do
            bcftools sort -O v "$file}" \
                > "${file/solution/sorted}"
            bcftools index "${file/solution/sorted}"
        done

        #bcftools "${data_dir}"/"${sample}"-*-MWER.sorted.vcf

    done

}

main
