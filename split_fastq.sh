#!/usr/bin/env bash

usage() {
  echo "Usage..."
}

# Default optional parameters
reads_per_file=10000000
threads=1
outdir="./"

while getopts ':n:t:o:' flag; do
  case "${flag}" in
    n)  reads_per_file="${OPTARG}" ;;
    t)  threads="${OPTARG}" ;;
    o)  outdir="${OPTARG}" ;;
    \?) usage
       exit 1 ;;
    :) echo "Invalid options: -"${OPTARG}" requires an argument" 1>&2
       exit 1 ;;
  esac
done

shift "$((OPTIND-1))"

if [[ ! -d "${outdir}" ]]; then
  >&2 echo "Error "${outdir}" is not a valid directory."
  exit 1
fi

for file in "${@}"; do
  if [[ ! -f "${file}" ]]; then
    >&2 echo "Error "${file}" is not a valid file."
    exit 1
  fi
done

split_fastq() {
  shopt -s extglob
  prefix="${outdir%*/}"/"${1##*/}"; prefix="${prefix%.f*q?(.gz)}".part
  zcat -f "${1}" \
    | split -l $(("${2}"*4)) --filter='gzip > "${FILE}".fq.gz' - "${prefix}"
}

export -f split_fastq

parallel -0 -j "${threads}" split_fastq ::: "${@}" ::: "${reads_per_file}" ::: "${outdir}"
