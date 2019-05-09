#!/usr/bin/env bash

usage() {
  echo "Usage: ..."
}

threads="1"
while getopts 'l:j:t:' flag; do
  case "${flag}" in
    l) ligation_site="${OPTARG}" ;;
    j) junction="${OPTARG}" ;;
	t) threads="${OPTARG}" ;;
    *) usage
       exit 1 ;;
  esac
done
shift "$((OPTIND-1))"

# Verify that mandatory parameters have been specified. 
if [ -z "${ligation_site}" ] || [ -z "${junction}" ]; then
    usage
	exit 1
fi

# Verify parameters are not files. Indicates error (i.e. truncate.sh -l GATCGATCG -j sample.fastq)
if [ -f "${ligation_site}" ] || [ -f "${junction}" ] || [ -f "${threads}" ]; then
    usage
	exit 1
fi

# Iterate through each file
for file in "$@"; do
	((i=i%threads)); ((i++==0)) && wait
	{	
  	echo $file
	if [ ! -f ${file} ]; then
    	echo "${file} not found."
	elif [[ ${file} == *.fastq ]]; then
		sed "/^>/! s/${ligation_site}.*/${junction}/" ${file} | gzip > ${file%.fastq}.trunc.fastq.gz
	elif [[ ${file} == *.fastq.gz ]]; then
		zcat ${file} | sed "/^>/! s/${ligation_site}.*/${junction}/" | gzip > ${file%.fastq.gz}.trunc.fastq.gz
	else
		echo "${file} is not fastq format (must end in .fastq or .fastq.gz)"
	fi
	} &
done; wait
