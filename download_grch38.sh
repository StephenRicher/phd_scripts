#!/usr/bin/env bash

#########################################################################
#### USER CONFIGURATION ####

# Specify genome build.
readonly genome_build=GRCh38
# Specficy genome URL
readonly genome_url=ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/\
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# Directory to store genomes.
readonly top_dir=~/x_db/DBuck/s_richer/genomes/

#########################################################################

# Directories
readonly genome_dir="${top_dir}"/"${build}"/
readonly qc_dir="${genome_dir}"/qc/
readonly idx_dir="${genome_dir}"/index/

# Genome index
readonly genome_idx="${idx_dir}"/"${build}"

# Files
readonly genome_path="${genome_dir}"/"${genome_url##*/}"

# Other
threads=6
main() {

    # Check top level directory exists.
    is_dir "${top_dir}" || fail "Error: "${top_dir}" is not a directory."

    # Create the necessary sub-directories.
    mkdir -p "${genome_dir}" "${qc_dir}" "${idx_dir}" || exit 1

    curl "${genome_url}" > "${genome_path}"

    bowtie2-build --threads "${threads}" --verbose \
        "${genome_path}" "${genome_idx}" \
        &> "${qc_dir}"/"${genome_build}".bt2_index.logfile

    bowtie2-inspect --summary "${genome_index}" \
        &> "${qc_dir}"/"${genome_build}".bt2_index_summary.txt
}

fail() {
    >&2 echo "${1}"
    exit "${2-1}"
}

main




