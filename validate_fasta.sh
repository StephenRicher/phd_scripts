#!/usr/bin/env bash

contains_whitespace() {
  # Returns true if file contains horizontal or vertical whitespace.
  if zcat -f "${1}" | grep -q -m 1 [[:space:]]; then
    return 0
  else
    return 1
  fi
}

invalid_fasta_start() {
  first_char=$(zcat -f "${1}" | grep -om 1 '^[[:blank:]]*.' | tr -d '[:blank:]')
  if [[ "${first_char}" == \> ]]; then
    return 1
  else
    return 0
  fi
}

fasta_contains_wrapping() {
  # Returns true if fasta sequences are wrapped (more than one sequence line per header).
  if zcat -f "${1}" | awk '/^[[:space:]]*>/ {n++; next} END {if(n*2 != NR) exit 1}'; then
    return 1
  else
    return 0
  fi
}

not_one_header_per_sequence() {
  # Extract first character of each header line, and line below, after removing blank lines and vertical whitespace.
  fasta_one_seq=( $(zcat -f "${1}" | sed '/^$/d' | tr -d '[:blank:]' | cut -c 1 | grep --no-group-separator -A 1 '^[[:blank:]]*>') )

  # Extract number of headers and sequences from array. Use null seperator in case header contains newline.
  num_headers=$(printf -- '%s\0' "${fasta_one_seq[@]}" | grep -cz '^>')
  num_sequences=$(printf -- '%s\0' "${fasta_one_seq[@]}" | grep -czv '^>')

  if [[ "${num_headers}" == "${num_sequences}" ]]; then
    return 1
  else
    return 0
  fi
}

contains_invalid_chars() {
  # Valid character set for nucleotide FASTA sequence
  valid_chars='^[ATCGURYKMSWBDHVN-]*$'

  # Isolate sequences, remove horiztonal whitespace and search for non valid characters.
  if zcat -f "${1}" | grep -v '^[[:blank:]]*>' | tr -d '[:blank:]' | grep -vicq -- "${valid_chars}"; then
    return 0
  else
    return 1
  fi
}


if contains_whitespace "${1}"; then
  >&2 echo "Error: Fasta file contains whitespace."
  exit 1
fi


if invalid_fasta_start "${1}"; then
  >&2 echo "Error: Fasta file does not begin with header"
  exit 1
fi

if fasta_contains_wrapping "${1}"; then
  >&2 echo "Error: Fasta file contains wrapping."
  exit 1
fi

if not_one_header_per_sequence "${1}"; then
  >&2 echo "Error: Number of fasta headers does not match number of sequences."
  exit 1
fi

if contains_invalid_chars "${1}"; then
  >&2 echo "Error: Fasta sequences contain invalid characters."
  exit 1
fi


