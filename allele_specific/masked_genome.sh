#!/usr/bin/env bash

shopt -s extglob

main() {
    local genome
    local vcf
    local dir="."
    local threads=1

    while getopts 'g:v:d:j:' flag; do
        case "${flag}" in
            g) genome="${OPTARG}" ;;
            v) vcf="${OPTARG}" ;;
            d) dir="${OPTARG%/}" ;;
            j) threads="${OPTARG}" ;;
            *) usage
            exit 1 ;;
        esac
    done

    if any_empty -n 2 "${genome}" "${vcf}"; then
       >&2 echo "Error: Missing mandatory arguments."
        usage
        exit 1
    fi

    if ! all_files "${genome}" "${vcf}"; then
        usage
        exit 1
    fi

    is_dir "${dir}" || exit 1

    # If genome is not indexed then index.
    if ! all_files "${genome}".fai; then
        samtools faidx "${genome}"
    fi

    ## Masked reference genome at SNPs and build bowtie2 index ##
    local genome_rmpath="${genome##*/}"
    local masked_genome="${dir}"/"${genome_rmpath%.fa?(.gz)}"_"${sample}"_masked.fa

    bedtools maskfasta \
        -fullHeader \
        -fi <(zcat -f "${genome}") \
        -bed "${vcf}" \
        -fo "${masked_genome}"

    local bt2_idx="${dir}"/GRCh38_"${sample}"
    bowtie2-build \
        --threads "${threads}" \
        "${masked_genome}" \
        "${dir}"/GRCh38_"${sample}"

    echo "${bt2_idx}"


usage() {
    >&2 echo usage...
}


main "${@}"
