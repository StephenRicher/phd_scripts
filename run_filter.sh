
digest_file="/home/stephen/x_db/DBuck/s_richer/stephen_test/genomes/GRCh38/digest/GRCh38_primary_assembly_Mbo1_digest.txt.gz"
filter_by_region="/home/stephen/h/phd/scripts2/hic_scripts/filter_by_region2.txt"
hicup_to_hicexplorer="/home/stephen/h/phd/scripts2/hic_scripts/hicup_to_hicexplorer.txt"
hicexplorer_normalize="/home/stephen/phd/scripts/hicexplorer_normalize.sh"
capture_regions="/home/stephen/h/phd/scripts2/hic_scripts/capture_regions.bed"
genome="/home/stephen/x_db/DBuck/s_richer/stephen_test/genomes/GRCh38/wgs/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
hicexplorer_normalize="/home/stephen/h/phd/scripts2/hic_scripts/hicexplorer_normalize.sh"
qc_dir="/home/stephen/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis/qc_logs/"
hictools="/home/stephen/phd/scripts/python3_modules/hictools.py"

for bam in *.coord_sorted.bam; do
  qualimap bamqc -bam "${bam}" \
                 --genome-gc-distr HUMAN \
                 --feature-file "${capture_regions}" \
                 --outside-stats \
                 -outdir "${qc_dir%*/}"/"${bam%%.*}"_bamQC \
                 --java-mem-size=4G
  echo -e ""${bam%%.*}"\t"${qc_dir%*/}"/"${bam%%.*}"_bamQC\t"${bam/_[12]*/}"" >> "${qc_dir%*/}"/multibamqc_config.txt
done
qualimap multi-bamqc --data "${qc_dir%*/}"/multibamqc_config.txt \
                     -outdir "${qc_dir%*/}"/multi_bamQC \
                     --java-mem-size=4G

for bam in *name_sorted.bam; do
  (
  sample="${bam%%.*}"; sample="${sample##*/}"
  #"${hictools}" process --digest "${digest_file}" --gunzip \
  #                       --output "${sample}".hic_info.bam "${bam}"
  # Sample 10% of reads and extract HiC stats.
  samtools view -h -s 42.05 "${sample}".hic_info.bam \
    | "${hictools}" extract --sample "${sample}" --gzip \
     >> all_samples2.hic_info_sampled.txt.gz
  )
done

# Remove all header from all_samples file and add 1 to the top of the file.
header=$(head -n 1 all_samples.hic_info_sampled.txt)
(printf "%s\n" "$header"; grep -vFxe "${header}" all_samples.hic_info_sampled.txt) > all_samples.hic_info_sampled.tmp.txt
mv all_samples.hic_info_sampled.tmp.txt all_samples.hic_info_sampled.txt

for bam in *hic_info.bam; do
  sample="${bam%%.*}"; sample="${sample##*/}"
  echo Processing "${bam}"
  samtools view -h "${bam}"\
    | "${filter_reads}" --min_inward 1000 --max_ditag 1000 \
    |  samtools sort > "${sample}".filtered-sort.bam
  samtools index "${sample}".filtered-sort.bam
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
  "${hicexplorer_normalize}" "${region}" "${chr}" "${start}" "${end}" "${hcx_dir}"/all_regions/"${region}" allele
done <${capture_regions}




