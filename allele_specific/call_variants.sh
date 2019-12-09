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
}


fail() {
    >&2 echo "${1}"
    exit "${2-1}"
}

usage() {
    >&2 echo usage...
}


main "${@}"

