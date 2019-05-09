#!/usr/bin/env bash

restriction_site="${1}"
overhang=$(echo "${restriction_site}" | grep -bo '\^' | cut -d ':' -f 1)
if [ -z "${overhang}" ]; then
  >&2 echo "Error: Restriction cut site not set. Use caret symbol to specify (e.g. ^GATC)."
fi
restriction_site=${restriction_site/^/}

fasta_contains_wrapping () {
  # Returns true if fasta sequences are wrapped (more than one sequence line per header).
  if zcat -f "${1}" | awk '/^[[:space:]]*>/ {n++; next} END {if(n*2 != NR) exit 1}'; then
    return 1
  else
    return 0
  fi
}

contains_whitespace () {
  # Returns true if file contains horizontal or vertical whitespace.
  if zcat -f "${1}" | grep -q -m 1 [[:space:]]; then
    return 0
  else
    return 1
  fi
}

if contains_whitespace "${2}"; then
  >&2 echo "Error: Fasta file contains whitespace."
  exit 1
fi

if fasta_contains_wrapping "${2}"; then
  >&2 echo "Error: Fasta file contains wrapping."
  exit 1
fi

while read line; do

  # If line is FASTA header.
  if [[ ${line} == \>* ]]; then
    ref_name=$(echo "${line}" | tr '[>\-_]' '\t' | cut -f 2)
    >&2 echo "Processing reference ${ref_name}."

  else
    ref_length=$(echo -n "${line}" | wc -c)
    declare -ia end_fragment=($(echo -n "${line}" | grep -ibo ${restriction_site} | cut -d ':' -f 1))
    # If no restriction sequence found then print whole sequence.
    if [ -z "${end_fragment}" ]; then
      printf "%s\t1\t%s\t1\n" "${ref_name}" "${ref_length}"

    else
      # If restriction sequence occurs at start, and there is no overhang, then remove.
      if [[ "${overhang}" == 0 ]] && [[ "${end_fragment[0]}" == 0 ]]; then
        unset end_fragment[0]
      fi
      declare -ia end_fragment=( "${end_fragment[@]/%/+"${overhang}"}" )
      declare -ia start_fragment=( "${end_fragment[@]/%/+1}" )
      start_fragment=("1" "${start_fragment[@]}")
      end_fragment+=( "${ref_length}" )
      ref_all=("${start_fragment[@]/*/"${ref_name}"}")
      # Index values
      declare -ia index=( "${!start_fragment[@]}" )
      declare -ia index=( "${index[@]/%/+1}" )
      paste <(printf "%s\n" "${ref_all[@]}") \
            <(printf "%s\n" "${start_fragment[@]}") \
            <(printf "%s\n" "${end_fragment[@]}") \
            <(printf "%s\n" "${index[@]}")
    fi
  fi
done < <(zcat -f "${2}")
