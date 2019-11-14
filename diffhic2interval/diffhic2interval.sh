#!/usr/bin/env bash

# Define directory of differential interactions.
dir=~/x_db/DBuck/s_richer/hic_01/data/diffhic/differential_interactions
# Define capture regions BED file.
capture_regions=/home/stephen/phd/scripts/capture_regions.bed

for comp in MCF7 HB2_CL4; do

    # For each comparison, generate seperate interval file for UP and DOWN.
    out_up="${dir}"/all_HB2_WT-vs-"${comp}"_di_up_interval.txt
    out_down="${dir}"/all_HB2_WT-vs-"${comp}"_di_down_interval.txt

    # Create header for output files and customise.
    for outfile in "${out_up}" "${out_down}"; do
        cp interval_header.txt "${outfile}"
        if [[ "${outfile}" == *"di_up"* ]]; then
            change="UP"
        else
            change="DOWN"
        fi
        # Modify header for each comparison and UP/DOWN logFC.
        sed -i -e "s/CHANGE/"${change}"/g" -e "s/COMP/"${comp}"/g" "${outfile}"
    done

    # Compute max logFC and scale scores by region.
    while IFS=$'\t' read -r chr start end region; do

        # Get maximum absolute logFC from UP/DOWN interactions.
        max=$(cat "${dir}"/"${region}"*"${comp}"*.arc \
              | cut -f 7 | tr -d '-' | sort | tail -n 1)

        cat "${dir}"/"${region}"*"${comp}"*.arc \
        | awk -v max="${max}" -v out_up="${out_up}" \
              -v out_down="${out_down}" -f interval_convert.awk

    done <"${capture_regions}"

done

