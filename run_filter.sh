filter="/home/stephen/phd/scripts/hic_filter/hic_filter.py"
digest_file="/media/stephen/Data/genomes/GRCh38/digest/GRCh38_primary_assembly_Mbo1_digest.txt"
filter_reads="/home/stephen/phd/scripts/hic_filter/filter_reads.py"
extract_hic_info="/home/stephen/phd/scripts/hic_filter/extract_hic_info.py"
filter_by_region="/home/stephen/phd/scripts/filter_by_region.txt"
hicup_to_hicexplorer="/home/stephen/phd/scripts/hicup_to_hicexplorer.txt"

capture_regions="/home/stephen/phd/scripts/capture_regions.bed"
genome="/media/stephen/Data/genomes/GRCh38/wgs/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

for bam in *.name_sorted.bam; do
  (
  sample="${bam%%.*}"; sample="${sample##*/}"
  #samtools view -h "${bam}" | "${filter}" -d "${digest_file}" - | samtools view -Sb > "${bam/name_sorted/hic_info}"
  # Sample 10% of reads and extract HiC stats.
  #samtools view -s 42.05 "${bam/name_sorted/hic_info}" | "${extract_hic_info}" --sample "${sample}" - >> all_samples.hic_info_sampled.txt
  ) &
done; wait

# Remove all header from all_samples file and add 1 to the top of the file.
#cat <(head -n 1 all_samples.hic_info_sampled.txt) <(grep -vP "sample\torientation" all_samples.hic_info_sampled.txt) > all_samples.hic_info_sampled.tmp.txt
#mv all_samples.hic_info_sampled.tmp.txt all_samples.hic_info_sampled.txt

for bam in *.name_sorted.bam; do
  echo Processing "${bam}"
  #samtools view -@ 12 -h "${bam/name_sorted/hic_info}" | "${filter_reads}" --min_inward 1000 --max_ditag 1000 - | samtools view -Sb > "${bam/name_sorted/filtered}"
  #samtools sort -@ 12 "${bam/name_sorted/filtered}" > "${bam/name_sorted/filtered-sort}"
  #samtools index -@ 12 "${bam/name_sorted/filtered-sort}"
done

#"${filter_by_region}" "${capture_regions}" "${genome}" 12 . *filtered-sort.bam

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
for bam in *filtered.bam; do
  ((count++))
  "${hicup_to_hicexplorer}" -f "${bam}" -c "${capture_regions}" -t 1 -o "${hcx_dir}" -d "${hcx_dir}"/"${genome_nopath%.*}"_"${enzyme}"_rest_site.bed

  # Combine all summary files. If first in loop then include header.
  if [ "${count}" -eq "1" ]; then
    cp  "${hcx_dir}"/"${bam%%.*}"-summary.txt "${hcx_dir}"/all_samples_summary.txt
  else
    tail -n +2 "${hcx_dir}"/"${bam%%.*}"-summary.txt >> "${hcx_dir}"/all_samples_summary.txt
  fi
done

