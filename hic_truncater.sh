#!/usr/bin/env bash

###################################################################################################################
# Truncates FASTQ sequences at ligation junction of restriction digest fragments. Accepts two positional arguments.
# First argument is restriction sequence with ^ to indicate cut site (e.g. Mbo1 = ^GATC).
# Second argument is a FASTQ file (either in GZ or uncompressed format).
# Truncated FASTQ file is written to stdout and summary stats are written to stderr.

invalid_fastq_start() {
  first_char=$(zcat -f "${1}" | grep -om 1 '^[[:blank:]]*.' | tr -d '[:blank:]')
  if [[ "${first_char}" == "@" ]]; then
    return 1  
  fi
}

if [[ $(echo "${1}" | grep -o '\^' | wc -l) != "1" ]]; then
  >&2 echo "Error: Restriction site "${1}" must contain one ^ to indicate cut site."
  exit 1
fi 

if [ ! -s "${2}" ]; then 
  >&2 echo "Error: "${2}" is not a valid file or is empty."
  exit 1 
fi

if invalid_fastq_start "${2}"; then
  >&2 echo "Error: "${2}" does not begin with FASTQ header (@ symbol)."
  exit 1
fi

# Extract restriction sequence (remove ^)
restriction_seq="${1/^/}"

# Calculate ligation sequence
cut_site1=$(echo "${1}" | grep -ob '\^' | cut -f 1 -d ':')
lig2=${restriction_seq:${cut_site1}}
cut_site2=$(("${#1}" -  ${cut_site1} - 1))
lig1=${restriction_seq:0:${cut_site2}}
ligation_seq=${lig1}${lig2}

# Remove all blank lines with sed before parsing with awk.
zcat -f "${2}" \
| sed '/^$/d' \
| awk -v f=${2} -v l=${ligation_seq} -v r=${replacement} '
  BEGIN {printf "file\ttotal_reads\ttruncated\tmean_truncated_sequence_length\n" > "/dev/stderr"}
  NR%4==2 {orig=length($0); gsub(l".*",r); len=length($0); if(orig!=len) {trunc+=1; trunc_len+=len}} 
  NR%4==0 {$0=substr($0,1,len)} {total+=1; print} 
  END {printf "%s\t%d\t%d\t%.2f\n",f,total/4,trunc,trunc_len/trunc > "/dev/stderr"}'
