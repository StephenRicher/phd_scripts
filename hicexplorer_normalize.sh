#!/usr/bin/env bash

# Activate extended globbing
shopt -s extglob

# Define script locations
deDoc="/home/stephen/h/phd/scripts2/hic_scripts/deDoc.jar"
ontad="/home/stephen/Downloads/OnTAD-master/src/OnTAD"
plot2dnSE="/home/stephen/h/phd/scripts2/hic_scripts/plot_twodnse.R"

region="${1}"
dir="${5}"
chr="${2}"
start="${3}"
end="${4}"

# Define plot region - set maximum plot region size
max_size=1500000
if (( $((end - start)) > "${max_size}" )); then
  mid=$(((end + start)/2))
  start_plot=$((mid-(max_size/2)))
  end_plot=$((mid+(max_size/2)))
  location="${chr}":"${start_plot}"-"${end_plot}"
else
  location="${chr}":"${start}"-"${end}"
fi

if [[ "${6}" == "allele" ]]; then
  samples=( "HB2_WT_G1" "HB2_WT_G2" "HB2_CL4_G1" "HB2_CL4_G2" "MCF7_G1" "MCF7_G2")
  bin_range=( $(seq 5000 1000 50000; seq 60000 5000 100000) )
  #tracks=("/home/stephen/x_am/RC-BB1219/stephen/projects/hic_analysis/hicexplorer_tracks/hic_track_genome_split.ini")
  tracks=("/home/stephen/x_am/RC-BB1219/stephen/projects/hic_analysis/hicexplorer_tracks/hic_track_genome_split_ontad.ini")
else
  samples=("HB2_WT" "HB2_CL4" "MCF7")
  bin_range=( $(seq 5000 1000 40000) rf)
  #tracks=("/home/stephen/x_am/RC-BB1219/stephen/projects/hic_analysis/hicexplorer_tracks/hic_track.ini" \
  #        "/home/stephen/x_am/RC-BB1219/stephen/projects/hic_analysis/hicexplorer_tracks/hic_track_replicate.ini" \
  #        "/home/stephen/x_am/RC-BB1219/stephen/projects/hic_analysis/hicexplorer_tracks/hic_track_replicate_ontad.ini")
  tracks=("/home/stephen/x_am/RC-BB1219/stephen/projects/hic_analysis/hicexplorer_tracks/hic_track_ontad.ini" \
          "/home/stephen/x_am/RC-BB1219/stephen/projects/hic_analysis/hicexplorer_tracks/hic_track_replicate_ontad.ini")
fi

# Check hic interaction density for each replicate of each sample. If any are 0 then remove from sample list
delete=()
for sample in "${samples[@]}"; do
  for rep in 1 2; do
    matrix="${dir}"/1000/"${sample}"-"${rep}"-"${region}"-1000.h5
    # Check if exists since empty directories will return null glob
    if [[ -f "${matrix}" ]]; then
      info=( $(hicInfo --matrices "${matrix}") )
      matrix_size=${info[13]/,/}
      matrix_sum=${info[19]/,/}
      matrix_density=$(echo "scale=3; ${matrix_sum} / (${matrix_size} * ${matrix_size})" | bc)
      if [[ "${matrix_density}" == "0" ]]; then
        >&2 echo "Removing "${sample}" from analysis due to insufficient interactions"
        delete+=("${sample}")
      fi
    fi
  done
done

for target in "${delete[@]}"; do
  for i in "${!samples[@]}"; do
    if [[ ${samples[i]} = "${target}" ]]; then
      unset 'samples[i]'
    fi
  done
done

echo "${samples[@]}"

# Of the retained samples, normalise the rf and 1000 bin size matrices seperately.
for bin in rf 1000; do
  # Find all matrices that match kept samples.
  matrices=()
  for sample in "${samples[@]}"; do
    while IFS=  read -r -d $'\0'; do
      matrices+=("$REPLY")
    done < <(find "${dir}"/"${bin}" -type f -name "${sample}*-"${bin}".h5" -print0)
  done
  matrices_norm=("${matrices[@]/.h5/-norm.h5}")
  hicNormalize --matrices "${matrices[@]}" --outFileName "${matrices_norm[@]}" --normalize smallest
done

# Note matrices_norm (used below) is defined for 1000 bin size matrices after for loop.
# Merge matrices at varying bin sizes and generate diagnostic plots to view bin read count distribution.
mkdir -p "${dir}"/correlation_plots "${dir}"/diagnostic_plots
rm "${dir}"/twodnSE.txt
for nbin in "${bin_range[@]}"; do

  # Skip if bin size is larger than the region
  if (( nbin > end - start )); then
    continue
  fi

  # Create new array so original is not overwritten
  samples2=( "${samples[@]}" )
  delete2=( "${delete[@]}" )

  # Merge each normalised and non-normalised matrix to the specified bin size.
  # Generate diagnostic plot for all matrices.
  mkdir -p "${dir}"/"${nbin}"
  
  for matrix in "${matrices_norm[@]}" "${matrices[@]}"; do
    matrix_rmpath="${matrix##*/}"; matrix_rmpath="${matrix_rmpath/'-1000'/-${nbin}}"
    if [[ "${nbin}" != "1000" ]] && [[ "${nbin}" != "rf" ]]; then
      hicMergeMatrixBins --matrix "${matrix}" --numBins $((${nbin}/1000)) --outFileName "${dir}"/"${nbin}"/"${matrix_rmpath}"
    fi
    hicCorrectMatrix diagnostic_plot --matrix "${dir}"/"${nbin}"/"${matrix_rmpath}" -o "${dir}"/diagnostic_plots/${matrix_rmpath%.h5}_diagnostic.png
    # Calculate onednSE for all matrices at all binsizes
    if [[ "${matrix}" != *"-norm.h5" ]] && [[ "${nbin}" != "rf" ]]; then
      hicConvertFormat --matrices "${dir}"/"${nbin}"/"${matrix_rmpath}" \
                       --outFileName "${dir}"/"${nbin}"/"${matrix_rmpath%.h5}".cool \
                       --inputFormat h5 --outputFormat cool 
      cooler dump --matrix "${dir}"/"${nbin}"/"${matrix_rmpath%.h5}".cool | awk '{$1=$1+1; $2=$2+1; print}' > "${dir}"/"${nbin}"/"${matrix_rmpath%.h5}".deDoc
      nrows=$(tail -n 1 "${dir}"/"${nbin}"/"${matrix_rmpath%.h5}".deDoc | cut -d ' ' -f 1)
      nbins=$((nrows * nrows))
      sed -i "1i\\${nbins}" "${dir}"/"${nbin}"/"${matrix_rmpath%.h5}".deDoc
      twodnSE=$(java -jar -Xmx8G "${deDoc}" "${dir}"/"${nbin}"/"${matrix_rmpath%.h5}".deDoc | tail -n 1 | rev | cut -d ' ' -f 1 | rev)
      echo -e "${matrix_rmpath%.h5}\t${region}\t${nbin}\t${twodnSE}" >> "${dir}"/twodnSE.txt
    fi

  done

  mkdir -p "${dir}"/"${nbin}"/tads
  # Merge replicate samples
  for sample in "${samples2[@]}"; do
    hicSumMatrices --matrices "${dir}"/"${nbin}"/"${sample}"*-norm.h5 --outFileName "${dir}"/"${nbin}"/"${sample}"-"${region}"-"${nbin}"-norm_sum.h5

    for matrix in "${dir}"/"${nbin}"/"${sample}"*-norm*(_sum).h5; do

      matrix_rmpath="${matrix##*/}"

      hicCorrectMatrix diagnostic_plot --matrix ${matrix} -o "${dir}"/diagnostic_plots/${matrix_rmpath%.h5}_diagnostic_temp.png \
        2> "${dir}"/diagnostic_plots/${matrix_rmpath%.h5}_diagnostic_log_temp.txt
      lower_threshold=$(grep "mad threshold" "${dir}"/diagnostic_plots/${matrix_rmpath%.h5}_diagnostic_log_temp.txt | awk '{print $NF}')
      hicCorrectMatrix correct --matrix "${matrix}" --correctionMethod ICE --iterNum 1000 \
                               --filterThreshold ${lower_threshold} 5 --outFileName "${matrix%.h5}"_iced.h5
   
      # Catch errors with ICE correct - remove sample if unsuccessful 
      if ! hicInfo --matrices "${matrix%.h5}"_iced.h5 2>/dev/null ; then
        unset 'samples2[sample]'
        delete2+=("${sample}")
        rm "${matrix%.h5}"_iced.h5
        continue
      fi

      rm "${dir}"/diagnostic_plots/${matrix_rmpath%.h5}_diagnostic_log_temp.txt 
      rm "${dir}"/diagnostic_plots/${matrix_rmpath%.h5}_diagnostic_temp.png
      if [[ "${nbin}" != "rf" ]]; then
        hicConvertFormat --matrices "${matrix%.h5}"_iced.h5 --inputFormat h5 --outputFormat homer \
                         --outFileName "${matrix%.h5}"_iced.homer.gz
        zcat "${matrix%.h5}"_iced.homer.gz | cut -f3- | tail -n+2 > "${matrix%.h5}"_iced.nn.txt
        "${ontad}" "${matrix%.h5}"_iced.nn.txt -maxsz $((1000000/nbin)) -o "${dir}"/"${nbin}"/tads/"${matrix_rmpath%.h5}"_iced -bedout chr"${chr}" "${nbin}"
        tail -n+2 "${dir}"/"${nbin}"/tads/"${matrix_rmpath%.h5}"_iced.bed \
          | awk -v start="${start}" -v OFS='\t' '
              {$2=$2+start; $3=$3+start; $7=$7+start; $8=$8+start; print $1, $2, $3, $1, $7, $8, 0}' \
          > "${dir}"/"${nbin}"/tads/"${matrix_rmpath%.h5}"_iced.links
        hicConvertFormat --matrices "${matrix%.h5}"_iced.h5 --inputFormat h5 --outputFormat cool \
                         --outFileName "${matrix%.h5}"_iced.cool --resolutions "${nbin}"
      fi
    done
  done

  echo "${samples2[@]}"
  if [ ${#samples2[@]} -eq 0 ]; then
    continue
  fi

  # Correlate all normalised ICED matrices
  hicCorrelate --matrices "${dir}"/"${nbin}"/*-norm_iced.h5 --outFileNameHeatmap "${dir}"/correlation_plots/${region}-"${nbin}"-heatmap.png \
               --outFileNameScatter "${dir}"/correlation_plots/${region}-"${nbin}"-scatter.png --threads 6 --method pearson


  # Iterate through replicate pairs and detect hierarchical TADs using hitad
  #if [[ "${nbin}" != "rf" ]]; then
  #  mkdir -p "${dir}"/"${nbin}"/tads
  #  for sample in "${samples[@]}"; do
  #    echo res:"${nbin}" > "${dir}"/"${nbin}"/tads/"${sample}"_hitad_metafile.txt
  #    r=1
  #    for replicate in "${dir}"/"${nbin}"/"${sample}"*-norm_iced.cool; do
  #      echo rep"${r}":"${replicate}" >> "${dir}"/"${nbin}"/tads/"${sample}"_hitad_metafile.txt
  #      ((r++))
  #    done
  #    hitad --output "${dir}"/"${nbin}"/tads/"${sample}"_hitad.txt \
  #          --datasets "${dir}"/"${nbin}"/tads/"${sample}"_hitad_metafile.txt \
  #          --logFile "${dir}"/"${nbin}"/tads/"${sample}"_hitad_log.txt
  #  done 
  #fi 

  mkdir -p "${dir}"/"${nbin}"/matrix_comparison/
  for sample1 in "${samples2[@]}"; do
    for sample2 in "${samples2[@]}"; do
      if [[ "${sample1}" == "${sample2}" ]]; then
        hicCompareMatrices --matrices "${dir}"/"${nbin}"/"${sample1}"*-norm_iced.h5 \
                       --outFileName "${dir}"/"${nbin}"/matrix_comparison/"${sample1}"-"${region}"_replicate-"${nbin}"_log2.h5 \
                       --operation log2ratio
        hicPlotMatrix --matrix "${dir}"/"${nbin}"/matrix_comparison/"${sample1}"-"${region}"_replicate-"${nbin}"_log2.h5 \
                  --outFileName "${dir}"/"${nbin}"/matrix_comparison/"${sample1}"-"${region}"_replicate-"${nbin}"_log2.png \
                  --colorMap bwr --region ${location} --dpi 300 --vMin -4 --vMax 4 \
                  --title "${sample1}"-"${region}"_replicate_log2
      else
        hicCompareMatrices --matrices "${dir}"/"${nbin}"/"${sample1}"*-norm_sum_iced.h5 "${dir}"/"${nbin}"/"${sample2}"*-norm_sum_iced.h5 \
                           --outFileName "${dir}"/"${nbin}"/matrix_comparison/"${sample1}"_vs_"${sample2}"-"${region}"-"${nbin}"_log2.h5 \
                           --operation log2ratio
        hicPlotMatrix --matrix "${dir}"/"${nbin}"/matrix_comparison/"${sample1}"_vs_"${sample2}"-"${region}"-"${nbin}"_log2.h5 \
                      --outFileName "${dir}"/"${nbin}"/matrix_comparison/"${sample1}"_vs_"${sample2}"-"${region}"-"${nbin}"_log2.png \
                      --colorMap bwr --region ${location} --dpi 300 --vMin -4 --vMax 4 \
                      --title ${sample1}_vs_${sample2}_log2
      fi
    done
  done

  # Perform loop detection on unmerged and summed HiC matrices and plot. If no loops detected then run without the loops
  mkdir "${dir}"/"${nbin}"/hic_plots

  # Define vMax for standardised colour scales in hicPlots - different vMax for summed and not summed matrices
  vMax_sum=$(hicInfo --matrices "${dir}"/"${nbin}"/*-norm_sum_iced.h5 | grep Maximum | cut -d ':' -f 2 | sort -r | head -n 1)
  vMax_not_sum=$(hicInfo --matrices "${dir}"/"${nbin}"/*-norm_iced.h5 | grep Maximum | cut -d ':' -f 2 | sort -r | head -n 1)


  for matrix in "${dir}"/"${nbin}"/*-norm_?(sum_)iced.h5; do
    matrix_rmpath="${matrix##*/}"
    loops_file="${dir}"/"${nbin}"/hic_plots/${matrix_rmpath/.h5}_loops.bedgraph
    hicDetectLoops --matrix ${matrix} -o "${loops_file}" \
                   --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 \
                   --pValuePreselection 0.05 --pValue 0.5 --peakInteractionsThreshold 20

    if [[ "${matrix}" == *"_sum_"* ]]; then
      vMax=${vMax_sum}
    else
      vMax=${vMax_not_sum}
    fi

    if [ -f "${loops_file}" ]; then
      hicPlotMatrix --matrix ${matrix} --outFileName "${dir}"/"${nbin}"/hic_plots/${matrix_rmpath/.h5}.png \
                    --colorMap PuRd --log1p --region ${location} --loops "${loops_file}" \
                    --title ${matrix_rmpath/.h5} --dpi 300 --vMin 1 --vMax ${vMax}
    else
      hicPlotMatrix --matrix ${matrix} --outFileName "${dir}"/"${nbin}"/hic_plots/${matrix_rmpath/.h5}.png \
                    --colorMap PuRd --log1p --region ${location} \
                    --title ${matrix_rmpath/.h5} --dpi 300 --vMin 1 --vMax ${vMax}
    fi
  done

  montage "${dir}"/"${nbin}"/hic_plots/*-norm_iced.png -geometry 1200x1200+1+1 -tile 2x "${dir}"/"${nbin}"/hic_plots/${region}_${nbin}_hic_plots.png
  montage "${dir}"/"${nbin}"/hic_plots/*-norm_sum_iced.png -geometry 1200x1200+1+1 -tile x1 "${dir}"/"${nbin}"/hic_plots/${region}-${nbin}-sum_hic_plots.png

  mkdir -p "${dir}"/"${nbin}"/tads
  mkdir -p "${dir}"/"${nbin}"/tad_plots/tracks

  for track in "${tracks[@]}"; do
      
    plot_track="${dir}"/"${nbin}"/tad_plots/tracks/"${track##*/}"-"${region}".ini

    sed "s/-region-/-${region}-/g" ${track} > "${plot_track}"
    sed -i "s/-binsize-/-${nbin}-/g" "${plot_track}"
    sed -i "s/-binsize-/-${nbin}-/g" "${plot_track}"
    full_dir="$(realpath "${dir}")"
    sed -i "s/directory/"${full_dir//\//\\/}"\/"${nbin}"/g" "${plot_track}"

    for target in "${delete2[@]}"; do
      sed -i "/Start "${target}"/,/End "${target}"/d" "${plot_track}"
    done

    if [[ "${track}" == *"replicate"* ]]; then
      plotname="${dir}"/"${nbin}"/tad_plots/${region}_replicate_${nbin}.png
    else
      plotname="${dir}"/"${nbin}"/tad_plots/${region}_${nbin}.png
    fi

    # Plot track first to identifiy max HiC value across all HiC matrices.
    hicPlotTADs --tracks "${plot_track}" --region ${location} --outFileName "${plotname}" --dpi 300 2> "${dir}"/"${nbin}"/tad_plots/${track##*/}_plot_log_temp.txt

    # Extract max value from hicPlotTADs outut and insert to new temp track.ini file.
    max_hic_value=$(grep "max values for track" "${dir}"/"${nbin}"/tad_plots/${track##*/}_plot_log_temp.txt | awk '{print $NF}' | sort -nr | head -n 1)
    sed -i "s/#max_value = none/max_value = ${max_hic_value}/g" "${plot_track}"
        
    # Replot the graph.
    hicPlotTADs --tracks "${plot_track}" --region ${location} --outFileName "${plotname}" --title ${region} --dpi 300 
  done

  # Remove temporary files
  rm "${dir}"/"${nbin}"/tad_plots/${track##*/}_plot_log_temp.txt
done

/usr/bin/Rscript "${plot2dnSE}" "${dir}"/twodnSE.txt

