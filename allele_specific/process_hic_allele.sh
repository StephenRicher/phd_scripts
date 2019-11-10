#!/usr/bin/env bash

# Set master project directory
project_dir="/home/stephen/x_db/DBuck/s_richer/hic_01"
# Set allele directory
allele_dir="${project_dir}"/allele_specific
# Set qc directory 
qc="${allele_dir}"/qc
mkdir -p "${qc}"
# Set data directory
data_dir="${project_dir}"/data
# Set capture regions
capture_regions="/home/stephen/phd/scripts/capture_regions.bed"
# Unmasked genome
unmasked_genome="/home/stephen/x_db/DBuck/s_richer/stephen_test/genomes/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
# Set genome digest (from unmasked genome
unmasked_genome_digest="/home/stephen/x_db/DBuck/s_richer/stephen_test/genomes/GRCh38/GRCh38_Mbo1-digest.txt.gz"
# Set masked genome index
masked_genome_index="${allele_dir}"/GRCh38_HB2_WT

# Define samples
samples=( HB2_WT-1 HB2_WT-2  )
allele_samples=( HB2_WT_G1-1 HB2_WT_G2-1 HB2_WT_G1-2 HB2_WT_G2-2 )

# Set VCF and snpSplit file names
vcf="HB2_WT_all_hapcompass_MWER_solution.txt.vcf"
snpsplit="HB2_WT_MWER_solution_snpsplit.txt"

# Download hapcompass solution from serve
#rsync -avzh sr467@balena.bath.ac.uk:/home/u/sr467/scratch/projects/HiC/allele_specific/"${vcf}" "${allele_dir}"/
#rsync -avzh sr467@balena.bath.ac.uk:/home/u/sr467/scratch/projects/HiC/allele_specific/"${snpsplit}" "${allele_dir}"/

# Mask reference genome and build custom bowtie2 index
~/phd/scripts/allele_specific/masked_genome.sh \
  -s HB2_WT \
  -g "${unmasked_genome}" \
  -v "${allele_dir}"/"${vcf}" \
  -t 6 \
  -o "${allele_dir}"/

alignment_stats="${qc}"/alignment_stats.txt
rm "${alignment_stats}"
printf 'sample\ttotal_mapped_pairs\tboth_pairs_unmapped\t' >> "${alignment_stats}"
printf 'r1_map_r2_unmapped\tr2_map_r1_unmapped\tduplicate_pairs\t' >> "${alignment_stats}"
printf 'quality_filtered_pairs\ttotal_pairs_retained\n' >> "${alignment_stats}"

# Map R1 and R2 reads to masked genome
for sample in "${samples[@]}"; do

    # File for intermediate logfile
    intermediate="${allele_dir}"/"${sample}".fixmate.bam

    hictools map \
      --index "${masked_genome_index}"  \
      --sample "${sample}" \
      --log "${qc}"/"${sample}".bowtie2.logfile \
      --intermediate "${intermediate}" \
      --sensitivity sensitive \
      "${data_dir}"/"${sample}"-R[14]-trim-trunc.fq.gz \
      2> "${qc}"/"${sample}".alignment_stats.txt \
    | hictools deduplicate \
        --log "${qc}"/"${sample}".dedup.logfile \
        2> "${qc}"/"${sample}".dedup_stats.txt \
    | hictools process \
        --digest "${unmasked_genome_digest}" \
        --log "${qc}"/"${sample}".process.logfile \
        > "${allele_dir}"/"${sample}".proc.bam

    hictools filter \
      --min_ditag 100 --max_ditag 1000 \
      --min_inward 1000 \
      --log "${qc}"/"${sample}".filter.logfile \
      --sample "${sample}" \
      "${allele_dir}"/"${sample}".proc.bam \
      > "${allele_dir}"/"${sample}".filt.bam \
      2>> "${qc}"/filter_statistics.tsv

    total_pairs=$(grep -m 1 'reads; of these:' "${qc}"/"${sample}".alignment_stats.txt | awk '{print $1}')
    both_pair_unmapped=$(samtools view -cf 12 "${intermediate}")
    r1_map_r2_unmap=$(samtools view -c -f 72 -F 4 "${intermediate}")
    r2_map_r1_unmap=$(samtools view -c -f 136 -F 4 "${intermediate}")
    unmapped_pairs=$(( both_pair_unmapped + r1_map_r2_ummap + r2_map_r1_unmap ))
    duplicate_pairs=$(grep -m 1 'DUPLICATE' "${qc}"/"${sample}".dedup_stats.txt | awk '{print $3}')
    total_kept_pairs=$(samtools view -cf 64 "${allele_dir}"/"${sample}".proc.bam)
    qual_filtered_pairs=$(( total_pairs - total_kept_pairs - duplicate_pairs - unmapped_pairs ))

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${sample}" "${total_pairs}" "${both_pair_unmapped}" \
      "${r1_map_r2_unmap}" "${r2_map_r1_unmap}" \
      "${duplicate_pairs}" "${qual_filtered_pairs}" \
      "${total_kept_pairs}" >> "${alignment_stats}"

    ~/phd/scripts/allele_specific/SNPsplit.pl \
      --snp_file "${snpsplit}" \
      --hic \
      "${allele_dir}"/"${sample}".filt.bam

    for a in 1 2; do
      samtools merge -n "${sample/-*/_G${a}-${sample/*-/}}".bam \
      "${sample}"*G"${a}"_[GU]["${a}"A].bam
    done

done

hcx_dir="${allele_dir}"/hicexplorer
mkdir -p "${hcx_dir}"
diffhic_dir="${allele_dir}"/diffhic
mkdir -p "${diffhic_dir}"

summary_file="${hcx_dir}"/all_samples_summary.txt
printf 'sample\tcapture_region\tvalid_hic_pairs\tregion_length\thic_pairs_per_kb\n' \
  > "${summary_file}"

# Remove custom genome if it exists to prevent appending to existing.
genome_no_path="${unmasked_genome##*/}"
custom_genome="${diffhic_dir}"/"${genome_no_path%.fa}".captured_regions.fa
rm "${custom_genome}"
# Create custom genome and rename FASTA header to region name.
while IFS=$'\t' read -r chr start end region; do
  samtools faidx "${unmasked_genome}" "${chr}":$((start+1))-"${end}" \
    | sed "1 s/^.*$/>${region}/" \
    >> "${custom_genome}"
done <"${capture_regions}"


for sample in "${allele_samples[@]}"; do

  samtools view -f 0x40 -b "${allele_dir}"/"${sample}".bam \
    > "${allele_dir}"/"${sample}".R1.bam

  samtools view -f 0x80 -b "${allele_dir}"/"${sample}".bam \
    > "${allele_dir}"/"${sample}".R2.bam

  # Extract SAM header and remove chromosome lines
  samtools view -H "${allele_dir}"/"${sample}".bam \
    | grep -v "@SQ" > "${diffhic_dir}"/"${sample}".custom_header.sam

  while IFS=$'\t' read -r chr start end region; do

    # For each region insert a modified chromosome header after line 1.
    region_length=$((end - start))
    sed -i "2i @SQ\tSN:${region}\tLN:${region_length}" \
      "${diffhic_dir}"/"${sample}".custom_header.sam

    sub_dir="${hcx_dir}"/all_regions/"${region}"/1000
    mkdir -p "${sub_dir}"

    hicBuildMatrix \
      --samFiles "${allele_dir}"/"${sample}".R1.bam \
                 "${allele_dir}"/"${sample}".R2.bam \
      --region ${chr}:$((${start}+1))-${end} \
      --binSize 1000 \
      --outFileName ${sub_dir}/${sample}-${region}-1000.h5 \
      --outBam ${sub_dir}/${sample}-${region}.bam \
      --QCfolder ${sub_dir}/${sample}-${region}-1000_QC \
      --skipDuplicationCheck \
      --threads 6

    # Correct SAM header to remove unused chromosomes.
    cat <(samtools view -H ${sub_dir}/${sample}-${region}.bam \
            | awk -v chr="@SQ\tSN:${chr}\t" '$0 ~ chr || /@HD/ || /@PG/') \
         <(samtools view -S ${sub_dir}/${sample}-${region}.bam) \
      | samtools view -Sb > ${sub_dir}/${sample}-${region}.tmp.bam
    mv ${sub_dir}/${sample}-${region}.tmp.bam ${sub_dir}/${sample}-${region}.bam

    # Generate modified region BAM file
    samtools view ${sub_dir}/${sample}-${region}.bam \
      | awk -v OFS='\t' -v chr="${chr}" -v len="${region_length}" -v start="${start}" -v region="${region}" '
          {$3=region; $4=$4-start+1; $8=$8-start+1} {print}' \
      >> "${diffhic_dir}"/"${sample}".captured.sam

    read_count=$(($(samtools view -c ${sub_dir}/${sample}-${region}.bam) / 2))
    hic_pairs_per_kb=$((${read_count} / (${region_length} / 1000)))
    printf '%s\t%s\t%s\t%s\t%s\t\n' \
      "${sample}" "${region}" "${read_count}" \
      "${region_length}" "${hic_pairs_per_kb}" \
      >> "${summary_file}"

  done <${capture_regions}

  cat "${diffhic_dir}"/"${sample}".custom_header.sam \
      "${diffhic_dir}"/"${sample}".captured.sam \
    | samtools view -Sb > "${diffhic_dir}"/"${sample}".captured.bam \

  rm "${allele_dir}"/"${sample}".R[12].bam
  rm "${diffhic_dir}"/"${sample}".captured.sam

done

# Run hicexplorer script
while IFS=$'\t' read -r chr start end region; do
  for binsize in $(seq 5000 1000 20000); do
    /home/stephen/phd/scripts/hicexplorer_normalize_v2.sh \
      -r "${region}" \
      -c "${chr}" \
      -s "${start}" \
      -e "${end}" \
      -b "${binsize}" \
      -a -t 4 \
      -d "${hcx_dir}"/all_regions/"${region}" "${allele_samples[@]}"
  done
done <"${capture_regions}"

# Switch to py27 environment for hifive
eval "$(conda shell.bash hook)"
conda activate py27env

# Run QUASAR-QC on all regions at multiple bin ranges
# Define function to join sequence by char
function join_by { local IFS="$1"; shift; echo "$*"; }

for sample in "${samples[@]}"; do
  matrix="${diffhic_dir}"/"${sample}".captured.bam
  hifive fends --length <(samtools view -H "${matrix}" | awk '$1 ~ /^@SQ/' | sed 's/:/ /g' | awk -v OFS='\t' '{print $3, $5}') \
               --binned 1000 -g GRCh38 "${sample}".fend
  hifive hic-data --bam <(samtools view -bf 0x40 "${matrix}") \
                        <(samtools view -bf 0x80 "${matrix}") \
                  --skip-duplicate-filtering \
                  "${sample}".fend \
                  "${sample}".hifive.mat
  hifive hic-project "${sample}".hifive.mat \
                     "${sample}".hifive.project
  while IFS=$'\t' read -r chr start end region; do
    hifive quasar -p "${sample}".hifive.project \
                  -r $(join_by , $(seq 1000 1000 50000)) \
                  -c "${region}" -o /dev/stdout -d 0 "${region}"-"${sample}".quasar \
      | head -n 53 | tail -n +4 \
      | awk -v region="${region}" -v sample="${sample}" '{print sample, $0, region}' \
      >> "${qc}"/quasar_qc.txt
    rm "${region}"-"${sample}".quasar
  done <"${capture_regions}"
  rm "${sample}".fend "${sample}".hifive.mat "${sample}".hifive.project
done


