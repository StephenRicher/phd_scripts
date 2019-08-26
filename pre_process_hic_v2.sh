#!/usr/bin/env bash

shopt -s extglob
export TMPDIR=/home/stephen/x_db/DBuck/s_richer/tmp/

## Set project directory ##
capture_regions="/home/stephen/phd/scripts/capture_regions.bed"
project=/home/stephen/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis2
qc="${project}"/qc
mkdir -p "${qc}"

data_dir="${project}"/data
mkdir -p "${data_dir}"

samples=( HB2_CL4-1 HB2_CL4-2 HB2_WT-1 HB2_WT-2 MCF7-1 MCF7-2 )

## Download reference genome ##
build=GRCh38
genome_dir=/home/stephen/x_db/DBuck/s_richer/stephen_test/genomes/"${build}"
mkdir -p "${genome_dir}"
genome_ftp=ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/\
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
genome="${genome_dir}"/"${genome_ftp##*/}"
curl "${genome_ftp}" > "${genome}"

## Generate restriction digest ##
genome_digest="${genome_dir}"/"${build}"_Mbo1-digest.txt.gz
hictools digest --log "${qc}"/"${build}"_Mbo1-digest.logfile \
  --restriction ^GATC -zu "${genome}" > "${genome_digest}"

## Generate bowtie2 index and inspect ## 
genome_index_dir="${genome_dir}"/bt2_index
mkdir -p "${genome_index_dir}"
genome_index="${genome_index_dir}"/"${build}"
bowtie2-build --threads 6 --verbose "${genome}" "${genome_index}" \
  &> "${qc}"/"${build}".bt2_index.logfile
bowtie2-inspect --summary "${genome_index}" \
  &> "${qc}"/"${build}".bt2_index_summary.txt

## Run FastQC on raw data ##
fastqc --threads 12 --outdir "${qc}" "${data_dir}"/*R[14].fastq.gz

## Run FastQ Screen on raw data ## UNTESTED####
fastq_screen --aligner bowtie2 --threads 6 --outdir "${qc}" \
    "${data_dir}"/*R[14].fastq.gz

## Remove adapter contamination with cutadapt ##
forward_adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
reverse_adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
parallel -j 12 \
    "cutadapt -a "${forward_adapter}" -A "${reverse_adapter}" \
    --error-rate 0.1 --overlap 1 --gc-content 46 \
    --minimum-length 20 --quality-cutoff 20 \
    --output "${data_dir}"/{1}-R1-trim.fq.gz \
    --paired-output "${data_dir}"/{1}-R4-trim.fq.gz \
    "${data_dir}"/{1}-R1.fastq.gz \
    "${data_dir}"/{1}-R4.fastq.gz \
    > "${qc}"/{1}.cutadapt.txt" \
  ::: "${samples[@]}"

## Truncate seqeuences at restriction ligation junction ##
parallel -j 12 \
    "hictools truncate --restriction ^GATC -zu {1} \
    > {=1 s/trim/trim-trunc/ =} \
    2>> "${qc}"/hictools-truncate_summary.txt" \
  ::: "${data_dir}"/*trim.fq.gz

## Run FastQC on adapter trimmed and truncated data ##
fastqc --threads 12 --outdir "${qc}" "${data_dir}"/*trim-trunc.fq.gz

alignment_stats="${qc}"/alignment_stats.txt
rm "${alignment_stats}"
printf 'sample\ttotal_mapped_pairs\tboth_pairs_unmapped\t' >> "${alignment_stats}"
printf 'r1_map_r2_unmapped\tr2_map_r1_unmapped\tduplicate_pairs\t' >> "${alignment_stats}"
printf 'quality_filtered_pairs\ttotal_pairs_retained\n' >> "${alignment_stats}"

# Map R1 and R2 reads
for sample in "${samples[@]}"; do

    # File for intermediate logfile
    intermediate="${data_dir}"/"${sample}".fixmate.bam

    hictools map \
      --index "${genome_index}"  \
      --sample "${sample}"
      --log "${qc}"/"${sample}".bowtie2.logfile \
      --intermediate "${intermediate}" \
      "${data_dir}"/"${sample}"-R[14]-trim-trunc.fq.gz \
      2> "${qc}"/"${sample}"_alignment_stats.txt \
    | hictools deduplicate \
        --log "${qc}"/"${sample}".dedup.logfile \
        2> "${qc}"/"${sample}".dedup_stats.txt \
    | hictools process \
        --digest "${genome_digest}" \
        --log "${qc}"/"${sample}".process.logfile \
        > "${data_dir}"/"${sample}".proc.bam

    total_pairs=$(grep -m 1 'reads; of these:' "${sample}".bowtie2_stats.txt | awk '{print $1}')
    both_pair_unmapped=$(samtools view -cf 12 "${intermediate}")
    r1_map_r2_unmap=$(samtools view -c -f 72 -F 4 "${intermediate}")
    r2_map_r1_unmap=$(samtools view -c -f 136 -F 4 "${intermediate}")
    unmapped_pairs=$(( both_pair_unmapped + r1_map_r2_ummap + r2_map_r1_unmap ))
    duplicate_pairs=$(grep -m 1 'DUPLICATE' "${sample}".dedup_stats.txt | awk '{print $3/2}')
    total_kept_pairs=$(samtools view -cf 64 "${sample}".proc.bam)
    qual_filtered_pairs=$(( total_pairs - total_kept_pairs - duplicate_pairs - unmapped_pairs ))

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${sample}" "${total_pairs}" "${both_pair_unmapped}" \
      "${r1_map_r2_unmap}" "${r2_map_r1_unmap}" \
      "${duplicate_pairs}" "${qual_filtered_pairs}" \
      "${total_kept_pairs}" >> "${alignment_stats}"

    hictools extract \
      --sample "${sample}" --gzip \
      --log "${qc}"/"${sample}".extract.logfile \
      <(samtools view -hs 42.05 "${data_dir}"/"${sample}".proc.bam) \
      >> "${qc}"/hic_filter_qc.txt.gz

done

# Remove any multiple headers in QC file
zcat "${qc}"/hic_filter_qc.txt.gz \
  | awk 'FNR==1 {header = $0; print} $0 != header' \
  | gzip > "${qc}"/hic_filter_qc.tmp.txt.gz
mv "${qc}"/hic_filter_qc.tmp.txt.gz "${qc}"/hic_filter_qc.txt.gz
~/phd/scripts/plot_filter.R "${qc}"/hic_filter_qc.txt.gz "${qc}"


hcx_dir="${data_dir}"/hicexplorer
mkdir -p "${hcx_dir}"
diffhic_dir="${data_dir}"/diffhic
mkdir -p "${diffhic_dir}"

summary_file="${hcx_dir}"/all_samples_summary.txt
printf 'sample\tcapture_region\tvalid_hic_pairs\tregion_length\thic_pairs_per_kb\n' \
  > "${summary_file}"

# Remove custom genome if it exists to prevent appending to existing.
genome_no_path="${genome##*/}"
custom_genome="${diffhic_dir}"/"${genome_no_path%.fa}".captured_regions.fa
rm "${custom_genome}"
# Create custom genome and rename FASTA header to region name.
while IFS=$'\t' read -r chr start end region; do
  samtools faidx "${genome}" "${chr}":$((start+1))-"${end}" \
    | sed "1 s/^.*$/>${region}/" \
    >> "${custom_genome}"
done <"${capture_regions}"

for sample in "${samples[@]}"; do

  hictools filter \
    --min_ditag 100 --max_ditag 1000 \
    --min_inward 1000 \
    --log "${qc}"/"${sample}".filter.logfile \
    "${data_dir}"/"${sample}".proc.bam \
    > "${data_dir}"/"${sample}".filt.bam

  samtools view -f 0x40 -b "${data_dir}"/"${sample}".filt.bam \
    > "${data_dir}"/"${sample}".R1.filt.bam

  samtools view -f 0x80 -b "${data_dir}"/"${sample}".filt.bam \
    > "${data_dir}"/"${sample}".R2.filt.bam

  # Extract SAM header and remove chromosome lines
  samtools view -H "${data_dir}"/"${sample}".filt.bam \
    | grep -v "@SQ" > "${diffhic_dir}"/"${sample}".custom_header.sam

  while IFS=$'\t' read -r chr start end region; do

    # For each region insert a modified chromosome header after line 1.
    region_length=$((end - start))
    sed -i "2i @SQ\tSN:${region}\tLN:${region_length}" \
      "${diffhic_dir}"/"${sample}".custom_header.sam

    sub_dir="${hcx_dir}"/all_regions/"${region}"/1000
    mkdir -p "${sub_dir}"

    hicBuildMatrix \
      --samFiles "${data_dir}"/"${sample}".R1.filt.bam \
                 "${data_dir}"/"${sample}".R2.filt.bam \
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

  done <${capture_regions}; wait

  cat "${diffhic_dir}"/"${sample}".custom_header.sam \
      "${diffhic_dir}"/"${sample}".captured.sam \
    | samtools view -Sb > "${diffhic_dir}"/"${sample}".captured.bam \

  rm "${data_dir}"/"${sample}".R[12].filt.bam
  rm "${diffhic_dir}"/"${sample}".captured.sam

done

while IFS=$'\t' read -r chr start end region; do
  for binsize in $(seq 1000 1000 10000); do
    /home/stephen/phd/scripts/hicexplorer_normalize_v2.sh -r "${region}" \
                                                          -c "${chr}" \
                                                          -s "${start}" \
                                                          -e "${end}" \
                                                          -b "${binsize}" \
                                                          -d "${hcx_dir}"/all_regions/"${region}" "${samples[@]}"
  done
done <"${capture_regions}"

