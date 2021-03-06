#!/usr/bin/env bash

main() {
    local vcf
    local start
    local end
    local sample
    local region
    local data_dir="."
    local threads=1

    while getopts 'v:c:s:e:n:r:d:j:' flag; do
        case "${flag}" in
            v) vcf="${OPTARG}" ;;
            c) chr="${OPTARG}" ;;
            s) start="${OPTARG}" ;;
            e) end="${OPTARG}" ;;
            n) sample="${OPTARG}" ;;
            r) region="${OPTARG}" ;;
            d) data_dir="${OPTARG}" ;;
            j) threads="${OPTARG}" ;;
            *) usage
            exit 1 ;;
        esac
    done

    if any_empty -n 6 "${vcf}" "${chr}" "${start}" \
                      "${end}" "${region}" "${sample}"; then
       >&2 echo "Error: Missing mandatory arguments."
        usage
        exit 1
    fi

    if ! all_files "${vcf}" "${vcf}".csi; then
        usage
        exit 1
    fi

    bcftools view \
            --output-type v \
            --regions "${chr}":$((${start}+1))-"${end}" \
            "${vcf}" \
        > "${data_dir}"/"${sample}"-"${region}".vcf
    java -Xmx100g -jar "${hapcompass}" \
        --bam "${data_dir}"/"${sample}".replicate_merge.bam \
        --vcf "${data_dir}"/"${sample}"-"${region}".vcf \
        --debug \
        --output "${data_dir}"/"${sample}"-"${region}"-hapcompass

    local MWER="${data_dir}"/"${sample}"-"${region}"-MWER_solution.txt
    local MWER_best=""${data_dir}"/"${sample}"-"${region}"-MWER_best.txt"

    # Get line number of header of top scoring section
    local top_score_line=$(grep -n BLOCK "${MWER}" \
                            | sort -k 6 -nr \
                            | cut -f 1 -d ':' \
                            | head -n 1)

    # Extract top scoring block from MWER solution
    awk -v n="${top_score_line}" '
            NR<n {next} NR==n {print;next} /^BLOCK/ {exit} {print}' \
        "${MWER}" \
        > "${MWER_best}"

    # Convert solution to VCF for genome masking
    java -jar "${hc2vcf}" "${MWER_best}" <(zcat "${vcf}") 2 true

    # Convert VCF to SNPsplit friendly output
    awk -v OFS=$'\t' '
            substr($10,1,3)=="0|1" {print $3, $1, $2, 1, $4"/"$5}
            substr($10,1,3)=="1|0" {print $3, $1, $2, 1, $5"/"$4}' \
            "${MWER_best}".vcf \
        > "${data_dir}"/"${sample}"-"${region}"-MWER_snpsplit.txt

    # Compress VCF
    bcftools view -O z "${MWER_best}".vcf > "${MWER_best}".vcf.gz
}


fail() {
    >&2 echo "${1}"
    exit "${2-1}"
}

usage() {
    >&2 echo usage...
}


main "${@}"
