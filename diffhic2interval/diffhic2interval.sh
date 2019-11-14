#!/usr/bin/env bash

dir=~/x_db/DBuck/s_richer/hic_01/data/diffhic/differential_interactions
capture_regions=/home/stephen/phd/scripts/capture_regions.bed

# Calculate score by region and by comparison.
for comp in MCF7 HB2_CL4; do

    # Create file with header to store all region for each comparison
    all_intervals="${dir}"/all_HB2_WT-vs-"${comp}"_di_all_interval.txt
    printf 'track type=interact name="HB2_WT-vs-'${comp}' DI CHANGE" ' > "${all_intervals}"
    printf 'description="Differential chromatin interactions: HB2_WT vs '${comp}' DI CHANGE" ' >> "${all_intervals}"
    printf 'interactUp=true useScore=on maxHeightPixels=200:100:50 visibility=full\n' >> "${all_intervals}"
    printf '#chrom\tchromStart\tchromEnd\tname\tscore\tvalue\texp\t' >> "${all_intervals}"
    printf 'color\tsourceChrom\tsourceStart\tsourceEnd\tsourceName\t' >> "${all_intervals}"
    printf 'sourceStrand\ttargetChrom\ttargetStart\ttargetEnd\t' >> "${all_intervals}"
    printf 'targetName\ttargetStrand\n' >> "${all_intervals}"

    while IFS=$'\t' read -r chr start end region; do

        out_all="${dir}"/"${region}"_HB2_WT-vs-"${comp}"_di_all.arc
        cat "${dir}"/"${region}"*"${comp}"*.arc > "${out_all}"

        if [ -s "${out_all}" ]; then

            interval="${dir}"/"${region}"_HB2_WT-vs-"${comp}"_di_all_interval.txt
            # Copy header to interval file
            head -n 2 "${all_intervals}" > "${interval}"
            max=$(cut -f 7 "${out_all}" | tr -d '-' | sort | tail -n 1)
            awk -v OFS='\t' -v max="${max}" '
                function abs(x){return ((x < 0.0) ? -x : x)}
                {score = (abs($7)*1000)/max}
                {if($7 < 0) colour = "#FF0000"; else colour = "#0000FF"}
                {print "chr"$1, $2, $6, ".", int(score), $7, ".", colour, "chr"$1, $2, $3, ".", ".", "chr"$4, $5, $6, ".", "."}' "${out_all}" >> "${interval}"

            # Append all record (except first two header lines to all_interval file
            tail -n +3 "${interval}" >> "${dir}"/all_HB2_WT-vs-"${comp}"_di_all_interval.txt
        else
            rm "${out_all}"
        fi
    done <"${capture_regions}"

    for change in "up" "down"; do
        all_intervals_change="${dir}"/all_HB2_WT-vs-"${comp}"_di_"${change}"_interval.txt
        if [ "${change}" == "up" ]; then
            sed 's/DI CHANGE/DI UP/g' "${all_intervals}" \
            | awk '(NR<3) || ($6 > 0)' > "${all_intervals/di_all/di_${change}}"
        else
            sed 's/DI CHANGE/DI DOWN/g' "${all_intervals}" \
            | awk '(NR<3) || ($6 < 0)' > "${all_intervals/di_all/di_${change}}"
        fi
    done

done

# Clean up intermediate files.
rm "${dir}"/*all_interval.txt


