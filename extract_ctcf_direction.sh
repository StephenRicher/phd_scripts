#!/usr/bin/env bash

print_usage() {
  echo "Usage: ..."
}

genome="${1}"
shift 1

if [ -z "${@}" ]; then
  print_usage
  exit 1
fi


ctcf_seq=$'>ctcf\nCCACNAGGTGGCAG'
ctcf_seq_revcomp=$'>ctcf\nCTGCCACCTNGTGG'

for file in "${@}"; do
    while IFS= read -r line; do
        if [[ "${line}" == "#"* ]]; then
            continue
        fi

        # Extract coordinates from input file
        ctcf=$(echo "${line}" | sed 's/chr//g' | awk '{print $2":"$3"-"$4}')

        declare -a alignments
        declare -a scores



        # Run alignment for forward and reverse sequence and save output.
        for sequence in "${ctcf_seq}" "${ctcf_seq_revcomp}"; do
            align=$(needle \
                -asequence <(samtools faidx "${genome}" "${ctcf}") \
                -bsequence <(echo "${sequence}") \
                -gapopen 10 -gapextend 0.5 \
                -outfile /dev/stdout 2> /dev/null)
            alignments+=("${align}")
            scores+=($(echo "${align}" | grep Score |  tr ' ' '\n' | tail -n 1))
        done

        # Skip if empty. Usually due to invalid chr name when running faidx.
        if [ -z "${alignments}" ]; then
            continue
        fi

        # Calculate alignment score difference between forward and reverse.
        score_diff=$(echo "${scores[0]} - ${scores[1]}"  | bc | tr -d '-')

        # Skip if score difference is less than 5 (directionality uncertain).
        if (( $(echo "${score_diff} < 5" | bc -l) )); then
            unset scores
            unset alignments
            continue
        elif (( $(echo "${scores[0]} > ${scores[1]}" | bc -l) )); then
            direction='+'
            score="${scores[0]}"
            alignment="${alignments[0]}"
        else
            direction='-'
            score="${scores[1]}"
            alignment="${alignments[1]}"
        fi

        # Extract only CTCF alignment line from 'needle' output.
        # Print each sequence character on its own line.
        # Get positions (line number) of sequence characters (not '-')
        # Extract start and end by sorting.
        ctcf_positions=$(echo "${alignment}" | grep 'ctcf ' \
        | awk -v ORS='' '{print $3}' \
        | sed 's/./\0\n/g' \
        | grep -nv '-' \
        | cut -d ':' -f 1)
        ctcf_start=$(echo "${ctcf_positions}" | sort -n | head -n 1)
        ctcf_end=$(echo "${ctcf_positions}" | sort -rn | head -n 1)

        # Extract genome start coordinate from 'asequence' of alignment output.
        genome_positions=$(echo "${alignment}" | grep '# 1:' | awk '{print $3}')
        genome_start="${genome_positions%-*}"

        # Calculate genome start/end coordinate of CTCF binding site
        abs_ctcf_start=$(("${genome_start}" + "${ctcf_start}" - 1))
        abs_ctcf_end=$(("${genome_start}" + "${ctcf_end}" - 1))

        unset scores
        unset alignments

        # Sort by column 1, then numeric sort by column 2
        echo "${line}" \
        | awk -v OFS='\t' -v direction="${direction}" \
              -v start="${abs_ctcf_start}" -v end="${abs_ctcf_end}" '
               {print $2, start, end, ".", $6, direction}'

    done < "${file}" | sort -k 1,1 -k 2,2n
done



