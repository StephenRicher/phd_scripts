#!/usr/bin/env bash

# Activate extended globbing
shopt -s extglob

while getopts 'r:c:s:e:b:d:at:' flag; do
  case "${flag}" in
    r) region="${OPTARG}" ;;
    c) chr="${OPTARG%/}" ;;
    s) start="${OPTARG}" ;;
    e) end="${OPTARG}" ;;
    b) binsize="${OPTARG}" ;;
    d) dir="${OPTARG}" ;;
    a) allele=true ;;
    t) threads="${OPTARG}" ;;
    *) print_usage 
       exit 1 ;;
  esac
done
shift "$((OPTIND-1))"

export NUMEXPR_MAX_THREADS="${threads}"

samples=("${@}")

# Get group name (i.e. without replicate) from sample name
groups=($(printf "%s\n" "${samples[@]/-*/}" | sort -u | tr '\n' ' '))

if [ "${allele}" = true ]; then
  tracks=("/home/stephen/phd/scripts/pyGenomeTracks_configs/hicexplorer_G1vsG2_ontad.ini" \
          "/home/stephen/phd/scripts/pyGenomeTracks_configs/hicexplorer_G1vsG2_ontad_sum.ini")
else
  tracks=("/home/stephen/phd/scripts/pyGenomeTracks_configs/hicexplorer_WTvsCL4vsMCF7_ontad.ini" \
          "/home/stephen/phd/scripts/pyGenomeTracks_configs/hicexplorer_WTvsCL4vsMCF7_ontad_sum.ini")
fi

# Define plot region - set maximum plot region size
max_size=2000000
if (( $((end - start)) > "${max_size}" )); then
  mid=$(((end + start)/2))
  start_plot=$((mid-(max_size/2)))
  end_plot=$((mid+(max_size/2)))
  plot_range="${chr}":"${start_plot}"-"${end_plot}"
  depth=$((max_size/2))
else
  plot_range="${chr}":"${start}"-"${end}"
  depth=$(((end - start)/2))
fi



# Check hic interaction density for each replicate of each sample. If any are 0 then remove from sample list
delete=()
for sample in "${samples[@]}"; do
  matrix="${dir}"/1000/"${sample}"-"${region}"-1000.h5
  # Check if exists (and not empty)
  if [[ -s "${matrix}" ]]; then
    info=( $(hicInfo --matrices "${matrix}") )
    matrix_size=${info[13]/,/}
    matrix_sum=${info[19]/,/}
    matrix_density=$(echo "scale=3; ${matrix_sum} / (${matrix_size} * ${matrix_size})" | bc)
    if [[ "${matrix_density}" == "0" ]]; then
      >&2 echo "Removing "${sample}" from analysis due to insufficient interactions"
      delete+=("${sample}")
    fi
  # If file exists but empty
  elif [[ -f "${matrix}" ]]; then
    delete+=("${sample}")
  fi
done

# Delete sample group (i.e. HB2_WT)
for target in "${delete[@]/-*/}"; do
  for i in "${!groups[@]}"; do
    if [[ ${groups[i]} = "${target}" ]]; then
      unset 'groups[i]'
    fi
  done
done

# Add retained matrix file locations to array and normalise.
matrices=()
for group in "${groups[@]}"; do
  while IFS=  read -r -d $'\0'; do
    matrices+=("$REPLY")
  done < <(find "${dir}"/1000 -type f -name "${group}*-1000.h5" -print0)
done
matrices_norm=("${matrices[@]/.h5/-norm.h5}")
hicNormalize --matrices "${matrices[@]}" \
             --outFileName "${matrices_norm[@]}" \
             --normalize smallest

mkdir -p "${dir}"/correlation_plots "${dir}"/diagnostic_plots "${dir}"/counts_vs_dist

# Skip if bin size is larger than the region
if (( binsize > end - start )); then
  exit 1
fi

# Create new array so original is not overwritten
#groups2=( "${groups[@]}" )
#delete2=( "${delete[@]}" )

mkdir -p "${dir}"/"${binsize}"

# Merge raw and normalised matrices into a variety of bin sizes
for matrix in "${matrices_norm[@]}" "${matrices[@]}"; do
  matrix_rmpath="${matrix##*/}" 
  matrix_rmpath="${matrix_rmpath/'-1000'/-${binsize}}"
  if [[ "${binsize}" != "1000" ]]; then
    hicMergeMatrixBins --matrix "${matrix}" \
                       --numBins $((${binsize}/1000)) \
                       --outFileName "${dir}"/"${binsize}"/"${matrix_rmpath}"
  fi
done

mkdir -p "${dir}"/"${binsize}"/tads
for group in "${groups[@]}"; do

  # Reformat to (n+3)*n matrix for hicep
  for matrix in "${dir}"/"${binsize}"/"${group}"*"${binsize}".h5; do
    hicConvertFormat --matrices "${matrix}" \
                     --inputFormat h5 --outputFormat homer \
                     --outFileName /dev/stdout \
      | zcat | sed 's/-/\t/g' | tail -n +2 \
      | awk -v OFS='\t' -v bin="${binsize}" -v chr="${region}" '
          {$0=gensub(/\s*\S+/,"",4); $1=chr; $2=(NR-1)*bin; $3=$2+bin; print}' \
      > "${matrix%.h5}"_hicrep.tsv
  done

  hicSumMatrices --matrices "${dir}"/"${binsize}"/"${group}"*-norm.h5 \
                 --outFileName "${dir}"/"${binsize}"/"${group}"-"${region}"-"${binsize}"-norm_sum.h5

  for matrix in "${dir}"/"${binsize}"/"${group}"*-norm*(_sum).h5; do
    matrix_rmpath="${matrix##*/}"
    hicCorrectMatrix diagnostic_plot --matrix ${matrix} \
                                     -o "${dir}"/diagnostic_plots/${matrix_rmpath%.h5}_diagnostic.png \
      2> "${dir}"/diagnostic_plots/${matrix_rmpath%.h5}_diagnostic_log.txt
    lower_threshold=$(grep "mad threshold" "${dir}"/diagnostic_plots/"${matrix_rmpath%.h5}"_diagnostic_log.txt | awk '{print $NF}')
    if [ -z "${lower_threshold}" ]; then
      delete+=("${group}")
      continue
    fi
    hicCorrectMatrix correct --matrix "${matrix}" --correctionMethod ICE --iterNum 1000 \
                             --filterThreshold "${lower_threshold}" 5 --outFileName "${matrix%.h5}"_iced.h5

    # Catch errors with ICE correct - remove sample if unsuccessful
    if ! hicInfo --matrices "${matrix%.h5}"_iced.h5 2>/dev/null ; then
      delete+=("${group}")
      rm "${matrix%.h5}"_iced.h5
      continue
    fi

    hicConvertFormat --matrices "${matrix%.h5}"_iced.h5 \
                     --inputFormat h5 --outputFormat homer \
                     --outFileName "${matrix%.h5}"_iced.homer.gz

    # Reformat to n*n matrix for OnTAD
    zcat "${matrix%.h5}"_iced.homer.gz | cut -f3- | tail -n+2 > "${matrix%.h5}"_iced.nn.txt
    OnTAD "${matrix%.h5}"_iced.nn.txt \
          -maxsz $((1000000/binsize)) \
          -o "${dir}"/"${binsize}"/tads/"${matrix_rmpath%.h5}"_iced \
          -bedout chr"${chr}" "${binsize}"
    tail -n+2 "${dir}"/"${binsize}"/tads/"${matrix_rmpath%.h5}"_iced.bed \
      | awk -v start="${start}" -v OFS='\t' '
          {$2=$2+start; $3=$3+start; $7=$7+start; $8=$8+start; print $1, $2, $3, $1, $7, $8, 0}' \
      > "${dir}"/"${binsize}"/tads/"${matrix_rmpath%.h5}"_iced.links
  done
done

# Remove groups that could not be iced
for target in "${delete[@]/-*/}"; do
  for i in "${!groups[@]}"; do
    if [[ ${groups[i]} = "${target}" ]]; then
      unset 'groups[i]'
    fi
  done
done

# Exit if no valid groups remain
if [ ${#groups[@]} -eq 0 ]; then
  exit 1
fi

# Correlate all normalised ICED matrices
hicCorrelate --matrices "${dir}"/"${binsize}"/*-norm_iced.h5 \
             --outFileNameHeatmap "${dir}"/correlation_plots/${region}-"${binsize}"-heatmap.png \
             --outFileNameScatter "${dir}"/correlation_plots/${region}-"${binsize}"-scatter.png \
             --threads 6 --method pearson


for matrix in "${dir}"/"${binsize}"/*-norm_?(sum_)iced.h5; do

  matrix_rmpath="${matrix##*/}"

  # Delete any empty file (see touc below) otherwise hicFindTADs breaks
  tad_score_bedgraph="${dir}"/"${binsize}"/tads/"${matrix_rmpath%.h5}"_tad_score.bm
  rm "${tad_score_bedgraph}"

  hicFindTADs --matrix "${matrix}" \
            --minDepth $((binsize*3)) \
            --maxDepth $((binsize*10)) \
            --step "${binsize}" \
            --outPrefix "${dir}"/"${binsize}"/tads/"${matrix_rmpath%.h5}" \
            --correctForMultipleTesting fdr

  # If no tads found create empty file so hicPlotTads doesn't break
  if [ ! -f "${tad_score_bedgraph}" ]; then
    touch "${tad_score_bedgraph}"
  fi

done


mkdir -p "${dir}"/"${binsize}"/matrix_comparison/
for ((i=0; i < ${#groups[@]}; i++)); do
  group1="${groups[i]}"
  hicCompareMatrices --matrices "${dir}"/"${binsize}"/"${group1}"*-norm_iced.h5 \
                     --outFileName "${dir}"/"${binsize}"/matrix_comparison/"${group1}"-"${region}"_replicate-"${binsize}"_log2.h5 \
                     --operation log2ratio
  hicPlotMatrix --matrix "${dir}"/"${binsize}"/matrix_comparison/"${group1}"-"${region}"_replicate-"${binsize}"_log2.h5 \
                --outFileName "${dir}"/"${binsize}"/matrix_comparison/"${group1}"-"${region}"_replicate-"${binsize}"_log2.png \
                --colorMap bwr --region "${chr}":"${start}"-"${end}" --dpi 300 --vMin -3 --vMax 3 \
                --title "${group1}"-"${region}"_replicate_log2

  for ((j=0; j < ${#groups[@]}; j++)); do
    group2="${groups[j]}"
    if [ "${group1}" == "${group2}" ]; then
      continue
    fi
    hicCompareMatrices --matrices "${dir}"/"${binsize}"/"${group1}"*-norm_sum_iced.h5 "${dir}"/"${binsize}"/"${group2}"*-norm_sum_iced.h5 \
                       --outFileName "${dir}"/"${binsize}"/matrix_comparison/"${group1}"_vs_"${group2}"-"${region}"-"${binsize}"_log2.h5 \
                       --operation log2ratio
    hicPlotMatrix --matrix "${dir}"/"${binsize}"/matrix_comparison/"${group1}"_vs_"${group2}"-"${region}"-"${binsize}"_log2.h5 \
                  --outFileName "${dir}"/"${binsize}"/matrix_comparison/"${group1}"_vs_"${group2}"-"${region}"-"${binsize}"_log2.png \
                  --colorMap bwr --region "${chr}":"${start}"-"${end}" --dpi 300 --vMin -3 --vMax 3 \
                  --title ${group1}_vs_${group2}-"${region}"_log2

    # Subtract tad insulation score from matrices between samples
    ~/phd/scripts/subtract-tads.py "${dir}"/"${binsize}"/tads/"${group1}"-"${region}"-"${binsize}"*tad_score.bm \
                                   "${dir}"/"${binsize}"/tads/"${group2}"-"${region}"-"${binsize}"*tad_score.bm \
    > "${dir}"/"${binsize}"/tads/"${group1}"-minus-"${group2}"-"${region}"-"${binsize}"-tad_score.bm

  done
done

# Create array of normalised iced matrices, extract sample names and generate dist vs counts plot
norm_iced=( "${dir}"/"${binsize}"/*-norm_iced.h5 )
norm_iced_sample=( ${norm_iced[@]##*/} )
norm_iced_sample=( ${norm_iced_sample[@]/-${region}*} )
hicPlotDistVsCounts --matrices "${norm_iced[@]}" \
  --plotFile  "${dir}"/counts_vs_dist/"${region}"_"${binsize}"_counts_vs_dist.png \
  --labels "${norm_iced_sample[@]}" \
  --maxdepth $(((end-start)/2)) \
  --plotsize 5 4.2

# Perform loop detection on unmerged and summed HiC matrices and plot. If no loops detected then run without the loops
mkdir "${dir}"/"${binsize}"/hic_plots
mkdir "${dir}"/"${binsize}"/loops
mkdir "${dir}"/"${binsize}"/hic_plots_obs_exp

# Define vMax for standardised colour scales in hicPlots - different vMax for summed and not summed matrices
vMax_sum=$(hicInfo --matrices "${dir}"/"${binsize}"/*-norm_sum_iced.h5 | grep Maximum | cut -d ':' -f 2 | sort -r | head -n 1)
vMax_not_sum=$(hicInfo --matrices "${dir}"/"${binsize}"/*-norm_iced.h5 | grep Maximum | cut -d ':' -f 2 | sort -r | head -n 1)

for matrix in "${dir}"/"${binsize}"/*-norm_?(sum_)iced.h5; do

  matrix_rmpath="${matrix##*/}"

  hicTransform -m "${matrix}" --method obs_exp -o "${dir}"/"${binsize}"/${matrix_rmpath/.h5}_obs_exp.h5

  loops_file="${dir}"/"${binsize}"/loops/${matrix_rmpath/.h5}_loops.bedgraph
  timeout 5m hicDetectLoops --matrix ${matrix} -o "${loops_file}" \
                 --minLoopDistance "${binsize}" --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 \
                 --pValuePreselection 0.05 --pValue 0.5 --peakInteractionsThreshold 20
  
  if [[ "$?" == 124 ]]; then
    >&2 echo "Error timeout on hicDetectLoops ${matrix}"
  fi

  if [[ "${matrix}" == *"_sum_"* ]]; then
    vMax=${vMax_sum}
  else
    vMax=${vMax_not_sum}
  fi

  if [ -f "${loops_file}" ]; then
    hicPlotMatrix --matrix ${matrix} --outFileName "${dir}"/"${binsize}"/hic_plots/${matrix_rmpath/.h5}.png \
                  --colorMap YlGn --log1p --region "${chr}":"${start}"-"${end}" --loops "${loops_file}" \
                  --title ${matrix_rmpath/.h5} --dpi 300 --vMin 1 --vMax ${vMax}
    hicPlotMatrix --matrix "${dir}"/"${binsize}"/${matrix_rmpath/.h5}_obs_exp.h5 \
                  --outFileName "${dir}"/"${binsize}"/hic_plots_obs_exp/${matrix_rmpath/.h5}_obs_exp.png \
                  --colorMap YlGn --region "${chr}":"${start}"-"${end}" --loops "${loops_file}" \
                  --title ${matrix_rmpath/.h5}_obs_exp --dpi 300 --vMin 0 --vMax 2
  else
    hicPlotMatrix --matrix ${matrix} --outFileName "${dir}"/"${binsize}"/hic_plots/${matrix_rmpath/.h5}.png \
                  --colorMap YlGn --log1p --region "${chr}":"${start}"-"${end}" \
                  --title ${matrix_rmpath/.h5} --dpi 300 --vMin 1 --vMax ${vMax}
    hicPlotMatrix --matrix "${dir}"/"${binsize}"/${matrix_rmpath/.h5}_obs_exp.h5 \
                  --outFileName "${dir}"/"${binsize}"/hic_plots_obs_exp/${matrix_rmpath/.h5}_obs_exp.png \
                  --colorMap YlGn --region "${chr}":"${start}"-"${end}" \
                  --title ${matrix_rmpath/.h5}_obs_exp --dpi 300 --vMin 0 --vMax 2
  fi
done

montage "${dir}"/"${binsize}"/hic_plots/*-norm_iced.png \
        -geometry +0+0 -tile 2x \
        "${dir}"/"${binsize}"/hic_plots/${region}_${binsize}_hic_plots.png
montage "${dir}"/"${binsize}"/hic_plots/*-norm_sum_iced.png \
        -geometry +0+0 -tile x1 \
        "${dir}"/"${binsize}"/hic_plots/${region}-${binsize}-sum_hic_plots.png
montage "${dir}"/"${binsize}"/hic_plots_obs_exp/*-norm_iced_obs_exp.png \
        -geometry +0+0 -tile 2x \
        "${dir}"/"${binsize}"/hic_plots_obs_exp/${region}_${binsize}_hic_plots_obs_exp.png
montage "${dir}"/"${binsize}"/hic_plots_obs_exp/*-norm_sum_iced_obs_exp.png \
        -geometry +0+0 -tile x1 \
        "${dir}"/"${binsize}"/hic_plots_obs_exp/${region}-${binsize}-sum_hic_plots_obs_exp.png

#for group in "${groups[@]}"; do
#  for type in "norm.h5" "norm_sum.h5"; do
#    echo res:"${binsize}" > "${dir}"/"${binsize}"/tads/"${group}"_hitad_metafile-"${type/.h5/}".txt
#    r=1
#    for replicate in "${dir}"/"${binsize}"/"${group}"*"${type}"; do
#      replicate_rmpath="${replicate##*/}"
#      out="${dir}"/"${binsize}"/"${replicate_rmpath%.h5}".cool
#      hicConvertFormat --matrices "${replicate}" \
#                       --outFileName "${out}" \
#                       --inputFormat h5 --outputFormat cool
#      nbins=$(($(tail -n 1 <(cooler dump "${out}") | cut -f 2) + 1))
#      chromsize=$((nbins*binsize))
#      # Load cooler and convert to normalised to integer
#      cooler load -f coo <(echo -e "${chr}\t${chromsize}"):"${binsize}" <(cooler dump "${out}" | awk -v OFS='\t' '$3=int($3)') "${out}".tmp.cool
#      mv "${out}".tmp.cool "${out}"
#      cooler balance --max-iters 1000 "${out}"
#      echo rep"${r}":"${out}" >> "${dir}"/"${binsize}"/tads/"${group}"_hitad_metafile-"${type/.h5/}".txt
#      ((r++))
#    done
#    hitad --output "${dir}"/"${binsize}"/tads/"${group}"-"${region}"-"${binsize}"-hitad-"${type/.h5/}".txt \
#          --datasets "${dir}"/"${binsize}"/tads/"${group}"_hitad_metafile-"${type/.h5/}".txt \
#          --logFile "${dir}"/"${binsize}"/tads/"${group}"_hitad_log-"${type/.h5/}".txt \
#          --maxsize 1000000
#    awk -v start="${start}" -v OFS='\t' '
#      {$2=$2+start; $3=$3+start; print $1, $2, $3, $1, $2, $3, 0}' \
#      "${dir}"/"${binsize}"/tads/"${group}"-"${region}"-"${binsize}"-hitad-"${type/.h5/}".txt \
#      > "${dir}"/"${binsize}"/tads/"${group}"-"${region}"-"${binsize}"-hitad-"${type/.h5/}".links
#  done
#done

mkdir -p "${dir}"/"${binsize}"/tad_plots/tracks

for transform in count obs_exp; do
  for track in "${tracks[@]}"; do

    track_rmpath="${track##*/}"
    plot_track="${dir}"/"${binsize}"/tad_plots/tracks/"${track_rmpath/.ini/}"-"${region}"-"${transform}".ini

    sed "s/-region-/-${region}-/g" ${track} > "${plot_track}"
    sed -i "s/depth = 1000000/depth = ${depth}/g" "${plot_track}"
    sed -i "s/-binsize-/-${binsize}-/g" "${plot_track}"
    sed -i "s/-binsize-/-${binsize}-/g" "${plot_track}"
    full_dir="$(realpath "${dir}")"
    sed -i "s/directory/"${full_dir//\//\\/}"\/"${binsize}"/g" "${plot_track}"

    for target in "${delete[@]}"; do
      sed -i "/Start "${target}"/,/End "${target}"/d" "${plot_track}"
    done

    if [[ "${track}" == *"_sum.ini" ]]; then
      if [[ "${track}" == *"hitad"* ]]; then
        plotname="${dir}"/"${binsize}"/tad_plots/${region}_hitad_sum_${binsize}_"${transform}".png
      else
        plotname="${dir}"/"${binsize}"/tad_plots/${region}_ontad_sum_${binsize}_"${transform}".png
      fi
    elif [[ "${track}" == *"hitad.ini" ]]; then
      plotname="${dir}"/"${binsize}"/tad_plots/${region}_hitad_${binsize}_"${transform}".png
    else
      plotname="${dir}"/"${binsize}"/tad_plots/${region}_ontad_${binsize}_${transform}.png
    fi

    if [[ "${transform}" == "obs_exp" ]]; then
      sed -i "s/\.h5/_${transform}\.h5/g" "${plot_track}"
      sed -i 's/min_value = 1/min_value = 0/g' "${plot_track}"
      sed -i 's/#colormap/colormap/g' "${plot_track}"
      sed -i '/transform = log1p/d' "${plot_track}"
      sed -i "s/#max_value = none/max_value = 2/g" "${plot_track}"
    fi

    # Plot track first to identifiy max HiC value across all HiC matrices.
    hicPlotTADs --tracks "${plot_track}" --region ${plot_range} --outFileName "${plotname}" --dpi 50 2> "${dir}"/"${binsize}"/tad_plots/${track##*/}_plot_log_temp.txt

    # Extract max value from hicPlotTADs outut and insert to new temp track.ini file.
    max_hic_value=$(grep "max values for track" "${dir}"/"${binsize}"/tad_plots/${track##*/}_plot_log_temp.txt | awk '{print $NF}' | sort -nr | head -n 1)
    sed -i "s/#max_value = none/max_value = ${max_hic_value}/g" "${plot_track}"

    # Replot the graph.
    hicPlotTADs --tracks "${plot_track}" --region ${plot_range} \
                --outFileName "${plotname}" \
                --title "${region}"' at '"$((binsize/1000))"'kb bin size' \
                --dpi 600

    rm "${dir}"/"${binsize}"/tad_plots/${track##*/}_plot_log_temp.txt
  done
done


