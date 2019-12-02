#!/usr/bin/env bash


main() {
    local input
    local genome
    local capture_regions
    local data_dir="."
    local threads=1

    while getopts 'i:g:d:c:j:' flag; do
        case "${flag}" in
            i) input="${OPTARG}" ;;
            g) genome="${OPTARG}" ;;
            c) capture_regions="${OPTARG}" ;;
            d) data_dir="${OPTARG%/}" ;;
            j) threads="${OPTARG}" ;;
            *) usage ;;
        esac
    done
    shift "$((OPTIND-1))"

    if any_empty -n 3 "${input}" "${genome}" "${capture_regions}"; then
       >&2 echo "Error: Missing mandatory arguments."
        usage
        exit 1
    fi

    if ! all_files "${input}" "${genome}" "${capture_regions}"; then
        usage
        exit 1
    fi

    hcx_dir="${data_dir}"/hicexplorer
    diffhic_dir="${data_dir}"/diffhic
    mkdir -p "${hcx_dir}" "${diffhic_dir}"

    local sample=$(get_sample "${input}") || exit 1

    # Split input into R1 and R2 reads
    parallel -j "${threads}" --xapply \
        "samtools view -f {1} -b "${input}" \
        > "${data_dir}"/"${sample}".{2}.filt.bam" \
        ::: 0x40 0x80 ::: R1 R2

    # Extract SAM header and remove chromosome lines
    samtools view -H "${input}" \
        | grep -v "@SQ" > "${diffhic_dir}"/"${sample}".custom_header.sam

    while IFS=$'\t' read -r chr start end region; do

        # For each region insert a modified chromosome header after line 1.
        region_length=$((end - start))
        sed -i "2i @SQ\tSN:${region}\tLN:${region_length}" \
            "${diffhic_dir}"/"${sample}".custom_header.sam

        sub_dir="${hcx_dir}"/all_regions/"${region}"/1000
        mkdir -p "${sub_dir}"

        local bam="${sub_dir}"/"${sample}"-"${region}".bam
        hicBuildMatrix \
                --samFiles "${data_dir}"/"${sample}".R1.filt.bam \
                           "${data_dir}"/"${sample}".R2.filt.bam \
                --region "${chr}":$(("${start}"+1))-"${end}" \
                --binSize 1000 \
                --outFileName "${sub_dir}"/"${sample}"-"${region}"-1000.h5 \
                --outBam "${bam}" \
                --QCfolder "${sub_dir}"/"${sample}"-"${region}"-1000_QC \
                --skipDuplicationCheck \
                --threads 6 \
            || continue

        remove_unused_chrs "${bam}" "${chr}"

        convert_to_hic "${bam}" "${sample}" "${region}"

        modify_bam "${bam}" "${region}" "${chr}" "${start}" \
            >> "${diffhic_dir}"/"${sample}".captured.sam

        count=$(($(samtools view -c ${sub_dir}/${sample}-${region}.bam) / 2))
        hic_pairs_per_kb=$((${count} / (${region_length} / 1000)))

        printf '%s\t%s\t%s\t%s\t%s\t\n%s\t%s\t%s\t%s\t%s\t\n' \
                "sample" "capture_region" "valid_hic_pairs" \
                "region_length" "hic_pairs_per_kb" \
               "${sample}" "${region}" "${count}" \
               "${region_length}" "${hic_pairs_per_kb}"

    done < "${capture_regions}"

    cat "${diffhic_dir}"/"${sample}".custom_header.sam \
            "${diffhic_dir}"/"${sample}".captured.sam \
        | samtools view -Sb \
        > "${diffhic_dir}"/"${sample}".captured.bam \

    rm "${data_dir}"/"${sample}".R[12].filt.bam \
       "${diffhic_dir}"/"${sample}".custom_header.sam \
       "${diffhic_dir}"/"${sample}".captured.sam
}

# Convert BAM file to .hic format HiC matrix for UCSC visualisation
convert_to_hic() {
    local bam="${1}"
    local sample="${2}"
    local region="${3}"
    local dir="${bam%/*}"

    # Convert BAM to TSV format for conversion to .hic format
    samtools view "${bam}" \
        | awk '
            /^@/ {next}
            {if(and($2,0x10)) s=0; else s=1}
            NR%2 {m=$5; printf "%s\t%s\t%s\t%s\t0\t", $1, s, $3, $4}
            NR%2==0 {printf "%s\t%s\t%s\t1\t%s\t%s\n", s, $3, $4, m, $5}' \
        > "${dir}"/"${sample}"-"${region}".pre.tsv

    # Convert TSV to .hic format
    juicer_tools pre \
        -r $(join_by , $(seq 1000 1000 20000)) \
        "${dir}"/"${sample}"-"${region}".pre.tsv \
        "${dir}"/"${sample}"-"${region}".hic \
        hg38


    # Merge tsv across replicates - juicer_tools doesn't work with fifo.
    local merged_pre="${dir}"/"${sample/-*/}"-sum.pre.tsv
    cat "${dir}"/"${sample/-*/}"*.pre.tsv > "${merged_pre}"
    # Create another .hic file merged across sample replicates
    juicer_tools pre \
        -r $(join_by , $(seq 1000 1000 20000)) \
        "${merged_pre}" \
        "${dir}"/"${sample/-*/}"-"${region}"-sum.hic \
        hg38
    rm "${merged_pre}"
}

# Remove chromosomes from SAM header except for the one passed to function.
remove_unused_chrs() {
    local bam="${1}"
    local chr="${2}"

    cat <(samtools view -H "${bam}" \
            | awk -v chr="@SQ\tSN:${chr}\t" '$0 ~ chr || /@HD/ || /@PG/') \
        <(samtools view -S "${bam}") \
        | samtools view -Sb > "${bam}".tmp.bam
    mv "${bam}".tmp.bam "${bam}"
}

# Adjust SAM to reference relative to the ref/start/end of a capture region.
modify_bam() {
    local bam="${1}"
    local region="${2}"
    local chr="${3}"
    local start="${4}"

    samtools view -h "${bam}" \
        | awk -v OFS='\t' -v chr="${chr}" \
              -v start="${start}" -v region="${region}" '
            {$3=region; $4=$4-start+1; $8=$8-start+1} {print}'
}

# Define function to join sequence by char
join_by() {
    local IFS="$1"
    shift
    echo "$*"
}


fail() {
    >&2 echo "${1}"
    exit "${2-1}"
}


usage() {
    >&2 echo usage...
}


main "${@}"


