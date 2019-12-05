#!/usr/bin/env bash


#########################################################################
#### USER CONFIGURATION ####

readonly project_name=hic-01
# Directory to store project folder.
readonly top_dir=/home/stephen/x_db/DBuck/s_richer/
# Genome
readonly build=GRCh38
readonly genome=/home/stephen/x_db/DBuck/s_richer/genomes/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
# Genome index
readonly grch38_idx=/home/stephen/x_db/DBuck/s_richer/genomes/GRCh38/bt2_index/GRCh38
# Define restriction enzyme cut sequence
readonly re_name=Mbo1
readonly re_seq='^GATC'
# Define capture regions
readonly capture_regions=/home/stephen/phd/scripts/capture_regions.bed
# Path to data URLs
readonly data_urls=data_urls.txt
# Define threads
threads=12


#########################################################################

# Directories
readonly project_dir="${top_dir}"/"${project_name}"/
readonly data_dir="${project_dir}"/data/
readonly qc_dir="${project_dir}"/qc
readonly paths_dir="${project_dir}"/paths/

# Files
readonly raw_fastqs="${paths_dir}"/raw_fastqs.txt
readonly trimmed_fastqs="${paths_dir}"/trimmed_fastqs.txt
readonly processed_bams="${paths_dir}"/hic-processed_bams.txt
readonly digest="${project_dir}"/"${build}"_"${re_name}"-digest.txt.gz

main() {

    # Check top level directory exists.
    all_dirs "${top_dir}" || fail

    # Check provided input files are not empty.
    all_files "${genome}" "${data_urls}" || fail

    # Create the necessary sub-directories.
    mkdir -p "${project_dir}" "${data_dir}" \
             "${qc_dir}" "${paths_dir}" \
        || fail

    # Create copy of data ULRS in project folder.
    cp "${data_urls}" "${paths_dir}"

    # Download data from URLs and write output paths to file.
    download_data -o "${data_dir}" -j "${threads}" "${data_urls}" \
        > "${raw_fastqs}" \
        || fail "Error: Non-zero return code in download_urls."

    fastqc --threads "${threads}" --outdir "${qc_dir}" $(< "${raw_fastqs}")

    fastq_screen --aligner bowtie2 --threads "${threads}" \
        --outdir "${qc_dir}" $(< "${raw_fastqs}")

    # Execute cutadapt
    local forward_adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    local reverse_adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    run_paired_cutadapt -i "${raw_fastqs}" -d "${data_dir}" \
        -q "${qc_dir}" -j "${threads}" -- \
        -a "${forward_adapter}" -A "${reverse_adapter}" \
        --error-rate 0.1 --overlap 1 --gc-content 46 \
        --minimum-length 20 --quality-cutoff 20 \
        > "${trimmed_fastqs}"

    # Generate genome digest file
    pyHiCTools digest --log "${qc_dir}"/"${build}"_"${re_name}"-digest.logfile \
                      --restriction "${re_seq}" -zu "${genome}" --verbose \
        > "${digest}"

    # Run pyHiCTools
    while read -r forward_reads; read -r reverse_reads; do
        run_pyHiCTools \
            -1 "${forward_reads}" \
            -2 "${reverse_reads}" \
            -x "${grch38_idx}" \
            -i "${digest}" \
            -s "${re_seq}" \
            -d "${data_dir}" \
            -q "${qc_dir}" \
            -j "${threads}" \
            -f
    done < "${trimmed_fastqs}" > "${processed_bams}"

    /home/stephen/phd/scripts/figures/plot_filter.R \
        <(cat "${data_dir}"/*.extracted-subsample.txt \
            | sed '2,${/sample\torientation/d;}') \
        "${qc_dir}"

    while read -r bam; do
        post_process_hic \
            -i "${bam}" \
            -g "${genome}" \
            -c "${capture_regions}" \
            -d "${data_dir}" \
            -j "${threads}"
    done < "${processed_bams}" \
        | sed '2,${/sample\tcapture_region/d;}' \
        > "${qc_dir}"/all_samples_summary.txt

    # Extract all samples names from processed bams file
    samples=( $(read_samples_from_file "${processed_bams}") )

    while IFS=$'\t' read -r chr start end region; do
        for binsize in $(seq 1000 1000 10000); do
            hicexplorer_normalize \
                -r "${region}" \
                -c "${chr}" \
                -s "${start}" \
                -e "${end}" \
                -b "${binsize}" \
                -t 4 \
                -d "${data_dir}"/hicexplorer/all_regions/"${region}" \
                "${samples[@]}"
        done
    done < "${capture_regions}"

}

read_samples_from_file() {
    local file="${1}"
    while read -r path; do
        get_sample "${path}"
    done < "${file}" | uniq
}

fail() {
    all_empty "${@}" || >&2 echo "${1}"
    exit "${2-1}"
}

main
