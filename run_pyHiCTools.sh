#!/usr/bin/env bash


main() {
    local foward
    local reverse
    local bt2_idx
    local digest
    local restriction_seq
    local data_dir="."
    local qc_dir="."
    local threads=1

    while getopts 'f:r:x:i:s:d:q:j:' flag; do
        case "${flag}" in
            f) forward="${OPTARG}" ;;
            r) reverse="${OPTARG}" ;;
            x) bt2_idx="${OPTARG}" ;;
            i) digest="${OPTARG}" ;;
            s) restriction_seq="${OPTARG}" ;;
            d) data_dir="${OPTARG%/}" ;;
            q) qc_dir="${OPTARG%/}" ;;
            j) threads="${OPTARG}" ;;
            *) usage ;;
        esac
    done
    shift "$((OPTIND-1))"

    if any_empty -n 5 "${forward}" "${reverse}" "${bt2_idx}" \
                      "${digest}" "${restriction_seq}"; then
       >&2 echo "Error: Missing mandatory arguments."
        usage
        exit 1
    fi

    if ! all_files "${forward}" "${reverse}" "${digest}"; then
        usage
        exit 1
    fi

    local sample=$(get_sample "${forward}") || exit 1

    intermediate="${data_dir}"/"${sample}".fixmate.bam
    hic_stats="${qc_dir}"/"${sample}".hic_stats.txt
    hic_processed="${data_dir}"/"${sample}".proc.bam
    hic_filtered="${data_dir}"/"${sample}".filt.bam
    hic_extract="${data_dir}"/"${sample}".extracted-subsample.txt
    truncation_summary="${qc_dir}"/"${sample}"-truncation_summary.txt
    forward_trunc=$(modify_path -d "${data_dir}" -a '-trunc' "${forward}")
    reverse_trunc=$(modify_path -d "${data_dir}" -a '-trunc' "${reverse}")

    # Check if any output files already exist.
    any_files "${intermediate}" "${hic_stats}" "${hic_processed}" \
              "${hic_filtered}" "${truncation_summary}" "${hic_extract}" \
              "${forward_trunc}" "${reverse_trunc}" \
        && exit 1

    # Run truncation in parallel - USE GNU PARALLEL WHEN THIS IS INSTALLED ON CLUSTER
    local i
    local input
    local output
    for fastq in forward reverse; do
        ((i=i%threads)); ((i++==0)) && wait
        (
        if [[ "${fastq}" == "forward ]]; then
            input="${forward}"
            output="${forward_trunc}"
        else
            input="${reverse}"
            output="${reverse_trunc}"
        pyHiCTools truncate --restriction "${restriction_seq}" -zu "${input}" \
            > "${output}" \
            2>> "${truncation_summary}"
        ) &
    done
    wait

    pyHiCTools map \
            --index "${bt2_idx}"  \
            --sample "${sample}" \
            --log "${qc_dir}"/"${sample}".bowtie2.logfile \
            --intermediate "${intermediate}" \
            -@ "${threads}" \
            "${forward_trunc}" "${reverse_trunc}" \
        2> "${qc_dir}"/"${sample}".alignment_stats.txt \
        | pyHiCTools deduplicate \
            --log "${qc_dir}"/"${sample}".dedup.logfile \
            -@ "${threads}" \
        2> "${qc_dir}"/"${sample}".dedup_stats.txt \
        | pyHiCTools process \
            --digest "${digest}" \
            --log "${qc_dir}"/"${sample}".process.logfile \
        > "${hic_processed}"

    total_pairs=$(grep -m 1 'reads; of these:' "${qc_dir}"/"${sample}".alignment_stats.txt | awk '{print $1}')
    both_pair_unmapped=$(samtools view -cf 12 "${intermediate}")
    r1_map_r2_unmap=$(samtools view -c -f 72 -F 4 "${intermediate}")
    r2_map_r1_unmap=$(samtools view -c -f 136 -F 4 "${intermediate}")
    unmapped_pairs=$(( both_pair_unmapped + r1_map_r2_ummap + r2_map_r1_unmap ))
    duplicate_pairs=$(grep -m 1 'DUPLICATE' "${qc_dir}"/"${sample}".dedup_stats.txt | awk '{print $3}')
    total_kept_pairs=$(samtools view -cf 64 "${hic_processed}")
    qual_filtered_pairs=$(( total_pairs - total_kept_pairs - duplicate_pairs - unmapped_pairs ))

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "sample" "total_mapped_pairs" "both_pairs_unmapped" \
        "r1_map_r2_unmapped" "r2_map_r1_unmapped" "duplicate_pairs" \
        "quality_filtered_pairs" "total_pairs_retained" \
        "${sample}" "${total_pairs}" "${both_pair_unmapped}" \
        "${r1_map_r2_unmap}" "${r2_map_r1_unmap}" \
        "${duplicate_pairs}" "${qual_filtered_pairs}" \
        "${total_kept_pairs}" > "${hic_stats}"

    samtools view -hs 42.05 "${hic_processed}" \
        | pyHiCTools extract \
            --sample "${sample}" \
            --log "${qc_dir}"/"${sample}".extract.logfile \
        >> "${hic_extract}"

    pyHiCTools filter \
        --min_ditag 100 --max_ditag 1000 \
        --min_inward 1000 \
        --log "${qc_dir}"/"${sample}".filter.logfile \
        --sample "${sample}" \
        "${hic_processed}" \
        > "${hic_filtered}" \
        2>> "${qc_dir}"/"${sample}"-filter_statistics.tsv

    echo "${hic_filtered}"
}


fail() {
    >&2 echo "${1}"
    exit "${2-1}"
}


usage() {
    >&2 echo usage...
}


main "${@}"
