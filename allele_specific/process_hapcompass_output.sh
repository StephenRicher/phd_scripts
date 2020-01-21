#!/usr/bin/env bash

main() {
    local data_dir="."

    while getopts 'v:d:' flag; do
        case "${flag}" in
            v) vcf="${OPTARG}" ;;
            d) data_dir="${OPTARG}" ;;
            *) usage
            exit 1 ;;
        esac
    done

    if any_empty -n 1 "${vcf}"; then
       >&2 echo "Error: Missing mandatory arguments."
        usage
        exit 1
    fi

    if ! all_files "${vcf}"; then
        usage
        exit 1
    fi

    local sample=$(get_sample "${vcf}") || exit 1

    MWER_solution_all=""${data_dir}"/"${sample}"_all_hapcompass_MWER_solution.txt"
    rm "${MWER_solution_all}"

    while IFS=  read -r -d $'\0' file; do
        # Get line number of header of top scoring section
        top_score_line=$(grep -n BLOCK "${file}" | sort -k 6 -nr | cut -f 1 -d ':' | head -n 1)
        # Print out only the top scoring block
        awk -v n="${top_score_line}" '
                NR<n {next} NR==n {print;next} /^BLOCK/ {exit} {print}' \
                "${file}" \
            >> "${MWER_solution_all}"

    done < <(find "${data_dir}" -type f -name ""${sample}"_*_hapcompass_MWER_solution.txt" -print0)

    # Convert solution to VCF for genome masking
    java -jar "${hc2vcf}" "${MWER_solution_all}" <(zcat "${vcf_out}") 2 true

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
