#!/usr/bin/env bash

shopt -s extglob

main() {
    local genome
    local capture_regions
    local input_bams
    local data_dir="."
    local threads=1

    while getopts 'g:c:d:j:' flag; do
        case "${flag}" in
            g) genome="${OPTARG}" ;;
            c) capture_regions="${OPTARG}" ;;
            d) data_dir="${OPTARG}" ;;
            j) threads="${OPTARG}" ;;
            *) usage
            exit 1 ;;
        esac
    done
    shift "$((OPTIND-1))"
    input_bams=("${@}")

    if any_empty -n 3 "${genome}" "${capture_regions}" "${input_bams}"; then
       >&2 echo "Error: Missing mandatory arguments."
        usage
        exit 1
    fi

    if ! all_files "${genome}" "${genome}".fai "${capture_regions}" \
                   "${input_bams[@]}"; then
        usage
        exit 1
    fi

    local sample=$(get_sample "${input_bams[1]}") || exit 1

    # Merge coordinate sorted replicates and index
    samtools merge -@ "${threads}" - "${input_bams[@]}" \
        | samtools sort -m 2G -@ "${threads}" \
        > "${data_dir}"/"${sample}".replicate_merge.bam

    samtools index -@ "${threads}" "${data_dir}"/"${sample}".replicate_merge.bam


    # Call variants from merged sample file and quality filter.
    vcf_out="${data_dir}"/"${sample}".sorted.vcf.gz

    bcftools mpileup \
            -q 15 --ignore-RG --count-orphans \
            --max-depth 100000 \
            --output-type u \
            -f "${genome}" \
            --regions-file "${capture_regions}" \
            "${data_dir}"/"${sample}".replicate_merge.bam \
        | bcftools call \
            --skip-variants indels \
            --multiallelic-caller \
            --variants-only \
            --output-type u \
        | bcftools view \
            -i '%QUAL>=20' \
            --output-type u \
        | bcftools sort \
            --output-file "${vcf_out}" \
            --output-type z

    bcftools index "${vcf_out}"

    # Subset VCF by capture region and perform haplotype phasing.
    while IFS=$'\t' read -r chr start end region; do
        bcftools view \
                --output-type v \
                --regions "${chr}":$(("${start}"+1))-"${end}" \
                "${vcf_out}" \
            > "${data_dir}"/"${sample}"_"${region}".vcf
        java -Xmx100g -jar "${hapcompass}" \
            --bam "${data_dir}"/"${sample}".replicate_merge.bam \
            --vcf "${data_dir}"/"${sample}"_"${region}".vcf \
            --debug \
            --output "${data_dir}"/"${sample}"_"${region}"_hapcompass
    done <"${capture_regions}"

    MWER_solution_all="${data_dir}/${sample}_all_hapcompass_MWER_solution.txt"
    rm "${MWER_solution_all}"
    while IFS=  read -r -d $'\0' file; do
        # Get line number of header of top scoring section
        top_score_line=$(grep -n BLOCK "${file}" | sort -k 6 -nr | cut -f 1 -d ':' | head -n 1)
        # Print out only the top scoring block
        awk -v n=${top_score_line} 'NR<n {next} NR==n {print;next} /^BLOCK/ {exit} {print}' ${file} >> ${MWER_solution_all}
    done < <(find "${data_dir}" -type f -name ""${sample}"_*_hapcompass_MWER_solution.txt" -print0)

    # Convert solution to VCF for genome masking
    java -jar "${hc2vcf}" "${MWER_solution_all}" \
        <(zcat "${vcf_out}") 2 true

    # Convert VCF to SNPsplit friendly output
    awk -v OFS=$'\t' 'substr($10,1,3)=="0|1" {print $3, $1, $2, 1, $4"/"$5} substr($10,1,3)=="1|0" {print $3, $1, $2, 1, $5"/"$4}' \
        "${MWER_solution_all}".vcf \
        > "${data_dir}"/"${sample}"_MWER_solution_snpsplit.txt

}


fail() {
    >&2 echo "${1}"
    exit "${2-1}"
}

usage() {
    >&2 echo usage...
}


main "${@}"

