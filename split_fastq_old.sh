#!/usr/bin/env bash

# Define cores
threads=16

# Specify number of reads to include per file.
reads=10000000

# Safely read file names into array using null character as seperator.
while IFS=  read -r -d $'\0' fastq; do
    zcat ${fastq} | split -l $((reads*4)) - ${fastq%.*}.part &
done < <(find $(pwd) -type f -name "*\.trunc\.fastq\.gz" -print0)
wait

# Gzip in parralel - 1 file per thread. Use -path because the search variable contains a path. (https://stackoverflow.com/questions/4341442/gzip-with-all-cores)
find $(pwd) -type f -path "*trunc\.fastq\.part*" -print0 | xargs -0 -n 1 -P ${threads} gzip

