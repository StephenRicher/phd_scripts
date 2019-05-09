#!/usr/bin/env bash

remove_fasta_wrapping () {
  # Remove all leading whitespace, replace whitespace in headers with '-' and remove whitespace in sequences.
  # Replace non-alphanumerics (except '>' and '-') with '_'.
  # Remove newline characters on lines not beginning with '>'.
  zcat -f "${1}" \
  | sed -e 's/^[[:space:]]*//' -e '/^>/s/[[:space:]]/-/g' -e 's/[[:blank:]]//' \
  | tr -c '[:alnum:][:space:].>-' '_' \
  | awk '!/^>/ {printf "%s", $0; n="\n"} /^>/ {print n$0; n=""} END {printf "%s", n}'
}

# Set input as first argument if defined, otherwise use standard input.
input="${1:-/dev/stdin}"

remove_fasta_wrapping "${input}"
