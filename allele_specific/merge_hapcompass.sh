#!/usr/bin/env bash

main() {
    local sample
    local samples=( 'HB2_CL4' 'HB2_WT' 'MCF7' )
    local data_dir="/home/u/sr467/scratch/projects/hic-01/allele_specific"

    for sample in "${samples[@]}"; do
        cat "${data_dir}"/"${sample}"*snpsplit* \
            > "${data_dir}"/"${sample}"_all_MWER_solution_snpsplit.txt

        bcftools concat --naive  \
                "${data_dir}"/"${sample}"*MWER*solution*vcf.gz \
            | bcftools view \
            > "${data_dir}"/"${sample}"-MWER_solution_all.vcf
    done
}

main
