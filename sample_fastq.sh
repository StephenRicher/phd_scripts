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

# If arg is '-' then shift to unset variable.
test "${1}" == "-" && shift
# If input unset or empty then set to /dev/stdin.
input="${1:-/dev/stdin}"

if [[ "${input}" != "/dev/stdin" ]]; then
  if [ -z "${input}" ]; then
    >&2 echo "Error: no file has been specified to sample."
    exit 1
  fi
  if [ ! -f "${input}" ]; then
    >&2 echo "Error: "${input}" is not a valid file."
    exit 1
  fi
fi

if [[ "${block_size}" == 1 ]]; then
  zcat -f "${input}" | shuf --random-source=<(get_seeded_random "${seed}") -n "${sample_size}"
else
  zcat -f "${input}" \
    | awk -v n="${block_size}" '
        ORS=NR%n?"\0":"\n"
        END {if(NR % n != 0) print "Error: Total lines not divisible by block size.\n" > "/dev/stderr"}' \
    | shuf --random-source=<(get_seeded_random "${seed}") -n "${sample_size}" \
    | tr "\0" "\n"
fi

