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

    while getopts '1:2:x:i:s:d:q:j:f' flag; do
        case "${flag}" in
            1) forward="${OPTARG}" ;;
            2) reverse="${OPTARG}" ;;
            x) bt2_idx="${OPTARG}" ;;
            i) digest="${OPTARG}" ;;
            s) restriction_seq="${OPTARG}" ;;
            d) data_dir="${OPTARG%/}" ;;
            q) qc_dir="${OPTARG%/}" ;;
            j) threads="${OPTARG}" ;;
            f) local keep_files="false" ;;
            *) usage ;;
        esac
    done
    shift "$((OPTIND-1))"


    any_empty -n 5 "${forward}" "${reverse}" "${bt2_idx}" \
                   "${digest}" "${restriction_seq}" \
        && fail "Error: Missing mandatory arguments."

    all_files "${forward}" "${reverse}" "${digest}" || fail

    all_dirs "${data_dir}" "${qc_dir}" || fail

    local sample=$(get_sample "${forward}") || fail

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
        && retain "${keep_files}" && fail

    # Run truncation in parallel - USE GNU PARALLEL WHEN THIS IS INSTALLED ON CLUSTER
    pyHiCTools truncate --restriction "${restriction_seq}" -zu "${forward}" \
            > "${forward_trunc}" 2>> "${truncation_summary}" & \
    pyHiCTools truncate --restriction "${restriction_seq}" -zu "${reverse}" \
            > "${reverse_trunc}" 2>> "${truncation_summary}"

    #export -f truncate
    #parallel -j "${threads}"  --xapply \
    #    truncate {1} {2} "${restriction_seq}" "${truncation_summary}" \
    #    ::: "${forward}" "${reverse}" \
    #    ::: "${forward_trunc}" "${reverse_trunc}"

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
    both_pair_unmapped=$(samtools view -c f 12 -@ "${threads}" "${intermediate}")
    r1_map_r2_unmap=$(samtools view -c -f 72 -F 4 -@ "${threads}" "${intermediate}")
    r2_map_r1_unmap=$(samtools view -c -f 136 -F 4 -@ "${threads}" "${intermediate}")
    unmapped_pairs=$(( both_pair_unmapped + r1_map_r2_ummap + r2_map_r1_unmap ))
    duplicate_pairs=$(grep -m 1 'DUPLICATE' "${qc_dir}"/"${sample}".dedup_stats.txt | awk '{print $3}')
    total_kept_pairs=$(samtools view -c f 64 -@ "${threads}" "${hic_processed}")
    qual_filtered_pairs=$(( total_pairs - total_kept_pairs - duplicate_pairs - unmapped_pairs ))

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "sample" "total_mapped_pairs" "both_pairs_unmapped" \
        "r1_map_r2_unmapped" "r2_map_r1_unmapped" "duplicate_pairs" \
        "quality_filtered_pairs" "total_pairs_retained" \
        "${sample}" "${total_pairs}" "${both_pair_unmapped}" \
        "${r1_map_r2_unmap}" "${r2_map_r1_unmap}" \
        "${duplicate_pairs}" "${qual_filtered_pairs}" \
        "${total_kept_pairs}" > "${hic_stats}"

    samtools view -h s 42.05 -@ "${threads}" "${hic_processed}" \
        | pyHiCTools extract \
            --sample "${sample}" \
            --log "${qc_dir}"/"${sample}".extract.logfile \
        >> "${hic_extract}"

    echo "${hic_processed}"
}

truncate() {
    local file="${1}"
    local output="${2}"
    local restriction="${3}"
    local summary="${4}"

    pyHiCTools truncate \
        --restriction "${restriction}" -zu "${file}" \
        >> "${output}" \
        2> "${summary}"
}


fail() {
    all_empty "${@}" || >&2 echo "${1}"
    usage
    exit "${2-1}"
}


usage() {
    >&2 echo usage...
}


main "${@}"
