#!/usr/bin/env bash

usage() {
  echo "Usage..."
}

# Default optional parameters
block_size=4
sample_size=1
seed="${RANDOM}"

while getopts ':b:n:s:' flag; do
  case "${flag}" in
    b)  block_size="${OPTARG}" ;;
    n)  sample_size="${OPTARG}" ;;
    s)  seed="${OPTARG}" ;;
    \?) usage
       exit 1 ;;
    :) echo "Invalid options: -"${OPTARG}" requires an argument" 1>&2
       exit 1 ;;
  esac
done

shift "$((OPTIND-1))"

get_seeded_random()
{
  seed="${1}"
  openssl enc -aes-256-ctr -pass pass:"${seed}" -nosalt \
    </dev/zero 2>/dev/null
}

if [ -z "${1}" ]; then
  >&2 echo "Error: no file has been specified to sample."
  exit 1
fi

if [ ! -f "${1}" ]; then
  >&2 echo "Error: "${1}" is not a valid file."
  exit 1
fi

if [[ "${block_size}" == 1 ]]; then
  zcat -f "${1}" | shuf --random-source=<(get_seeded_random "${seed}") -n "${sample_size}"
else
  zcat -f "${1}" \
    | awk -v n="${block_size}" '
        {printf("%s%s",$0,(NR%n==0)?"\n":"\0")}
        END {if(NR % n != 0) print "Error: Total lines not divisible by block size." > "/dev/stderr"}' \
    | shuf --random-source=<(get_seeded_random "${seed}") -n "${sample_size}" \
    | tr "\0" "\n"
fi

