filter="/home/stephen/h/phd/scripts2/hic_scripts/hic_filter/hic_filter.py"
digest_file="/home/stephen/x_db/DBuck/s_richer/stephen_test/genomes/GRCh38/digest/GRCh38_primary_assembly_Mbo1_digest.txt"
filter_reads="/home/stephen/h/phd/scripts2/hic_scripts/hic_filter/filter_reads.py"
extract_hic_info="/home/stephen/h/phd/scripts2/hic_scripts/hic_filter/extract_hic_info.py"
filter_by_region="/home/stephen/h/phd/scripts2/hic_scripts/filter_by_region.txt"
hicup_to_hicexplorer="/home/stephen/h/phd/scripts2/hic_scripts/hicup_to_hicexplorer.txt"
hicexplorer_normalize="/home/stephen/phd/scripts/hicexplorer_normalize.sh"
capture_regions="/home/stephen/h/phd/scripts2/hic_scripts/capture_regions.bed"
genome="/home/stephen/x_db/DBuck/s_richer/stephen_test/genomes/GRCh38/wgs/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
hicexplorer_normalize="/home/stephen/h/phd/scripts2/hic_scripts/hicexplorer_normalize.sh"

for bam in "${@}"; do
  (
  sample="${bam%%.*}"; sample="${sample##*/}"
  samtools view -h "${bam}" | "${filter}" -d "${digest_file}" - | samtools view -Sb > "${sample}".hic_info.bam
  # Sample 10% of reads and extract HiC stats.
  samtools view -s 42.05 "${sample}".hic_info.bam | "${extract_hic_info}" --sample "${sample}" - >> all_samples.hic_info_sampled.txt
  ) &
done; wait

# Remove all header from all_samples file and add 1 to the top of the file.
header=$(head -n 1 all_samples.hic_info_sampled.txt)
(printf "%s\n" "$header"; grep -vFxe "${header}" all_samples.hic_info_sampled.txt) > all_samples.hic_info_sampled.tmp.txt
mv all_samples.hic_info_sampled.tmp.txt all_samples.hic_info_sampled.txt

for bam in "${@}"; do
  sample="${bam%%.*}"; sample="${sample##*/}"
  echo Processing "${sample}".hic_info.bam
  samtools view -@ 6 -h "${sample}".hic_info.bam \
    | "${filter_reads}" --min_inward 1000 --max_ditag 1000 - \
    | samtools view -@ 6 -b > "${sample}".filtered.bam
  samtools sort -@ 3 "${sample}".filtered.bam > "${sample}".filtered-sort.bam
  samtools index -@ 3 "${sample}".filtered-sort.bam
done

"${filter_by_region}" "${capture_regions}" "${genome}" 12 . *filtered-sort.bam

# Specify restriction site and enzyme name.
enzyme="MboI"; restriction_site="GATC"

hcx_dir="./hicexplorer"
mkdir -p "${hcx_dir}"

genome_nopath="${genome##*/}"
# If digest file does not exist then create one.
if [ ! -f "${hcx_dir}"/"${genome_nopath%.*}"_"${enzyme}"_rest_site.bed ]; then
    findRestSite --fasta <(zcat -f "${genome}") \
                 --searchPattern "${restriction_site}" \
                 --outFile "${hcx_dir}"/"${genome_nopath%.*}"_"${enzyme}"_rest_site.bed
fi

count=0
for bam in *filtered-sort.bam; do
  ((count++))
  "${hicup_to_hicexplorer}" -f "${bam}" -c "${capture_regions}" -t 3 -o "${hcx_dir}" -d "${hcx_dir}"/"${genome_nopath%.*}"_"${enzyme}"_rest_site.bed

  # Combine all summary files. If first in loop then include header.
  if [ "${count}" -eq "1" ]; then
    cp  "${hcx_dir}"/"${bam%%.*}"-summary.txt "${hcx_dir}"/all_samples_summary.txt
  else
    tail -n +2 "${hcx_dir}"/"${bam%%.*}"-summary.txt >> "${hcx_dir}"/all_samples_summary.txt
  fi
done


while IFS=$'\t' read -r chr start end region; do
  "${hicexplorer_normalize}" "${region}" "${chr}" "${start}" "${end}" "${hcx_dir}"/all_regions/"${region}"
done <${capture_regions}




