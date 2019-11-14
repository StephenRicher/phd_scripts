#!/usr/bin/env bash

# Move to hic_analysis environment
eval "$(conda shell.bash hook)"
conda activate hic_analysis

shopt -s extglob
export TMPDIR=/home/stephen/x_db/DBuck/s_richer/tmp/

# Define function to join sequence by char
function join_by { local IFS="$1"; shift; echo "$*"; }

## Set project directory ##
capture_regions="/home/stephen/phd/scripts/capture_regions.bed"
project=/home/stephen/x_db/DBuck/s_richer/hic_01
qc="${project}"/qc
mkdir -p "${qc}"

data_dir="${project}"/data
mkdir -p "${data_dir}"

samples=( HB2_WT-1 HB2_WT-2 HB2_CL4-1 HB2_CL4-2 MCF7-1 MCF7-2 )

## Download reference genome ##
build=GRCh38
genome_dir=/home/stephen/x_db/DBuck/s_richer/stephen_test/genomes/"${build}"
mkdir -p "${genome_dir}"
genome_ftp=ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/\
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
genome="${genome_dir}"/"${genome_ftp##*/}"
genome="${genome%.gz}"
curl "${genome_ftp}" | zcat -f > "${genome}"

## Generate restriction digest ##
genome_digest="${genome_dir}"/"${build}"_Mbo1-digest.txt.gz
pyHiCTools digest --log "${qc}"/"${build}"_Mbo1-digest.logfile \
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
    "pyHiCTools truncate --restriction ^GATC -zu {1} \
    > {=1 s/trim/trim-trunc/ =} \
    2>> "${qc}"/pyHiCTools-truncate_summary.txt" \
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

    pyHiCTools map \
      --index "${genome_index}"  \
      --sample "${sample}" \
      --log "${qc}"/"${sample}".bowtie2.logfile \
      --intermediate "${intermediate}" \
      "${data_dir}"/"${sample}"-R[14]-trim-trunc.fq.gz \
      2> "${qc}"/"${sample}".alignment_stats.txt \
    | pyHiCTools deduplicate \
        --log "${qc}"/"${sample}".dedup.logfile \
        2> "${qc}"/"${sample}".dedup_stats.txt \
    | pyHiCTools process \
        --digest "${genome_digest}" \
        --log "${qc}"/"${sample}".process.logfile \
        > "${data_dir}"/"${sample}".proc.bam

    total_pairs=$(grep -m 1 'reads; of these:' "${qc}"/"${sample}".alignment_stats.txt | awk '{print $1}')
    both_pair_unmapped=$(samtools view -cf 12 "${intermediate}")
    r1_map_r2_unmap=$(samtools view -c -f 72 -F 4 "${intermediate}")
    r2_map_r1_unmap=$(samtools view -c -f 136 -F 4 "${intermediate}")
    unmapped_pairs=$(( both_pair_unmapped + r1_map_r2_ummap + r2_map_r1_unmap ))
    duplicate_pairs=$(grep -m 1 'DUPLICATE' "${qc}"/"${sample}".dedup_stats.txt | awk '{print $3}')
    total_kept_pairs=$(samtools view -cf 64 "${data_dir}"/"${sample}".proc.bam)
    qual_filtered_pairs=$(( total_pairs - total_kept_pairs - duplicate_pairs - unmapped_pairs ))

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${sample}" "${total_pairs}" "${both_pair_unmapped}" \
      "${r1_map_r2_unmap}" "${r2_map_r1_unmap}" \
      "${duplicate_pairs}" "${qual_filtered_pairs}" \
      "${total_kept_pairs}" >> "${alignment_stats}"

    pyHiCTools extract \
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

  pyHiCTools filter \
    --min_ditag 100 --max_ditag 1000 \
    --min_inward 1000 \
    --log "${qc}"/"${sample}".filter.logfile \
    --sample "${sample}" \
    "${data_dir}"/"${sample}".proc.bam \
  > "${data_dir}"/"${sample}".filt.bam \
  2>> "${qc}"/filter_statistics.tsv

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
    mv ${sub_dir}/${sample}-${region}.tmp.bam \
       ${sub_dir}/${sample}-${region}.bam

    # Convert BAM to TSV format for conversion to .hic format
    samtools view ${sub_dir}/${sample}-${region}.bam \
    | awk '
        /^@/ {next}
        {if(and($2,0x10)) s=0; else s=1}
        NR%2 {m=$5; printf "%s\t%s\t%s\t%s\t0\t", $1, s, $3, $4}
        NR%2==0 {printf "%s\t%s\t%s\t1\t%s\t%s\n", s, $3, $4, m, $5}' \
    > ${sub_dir}/${sample}-${region}.pre.tsv

    # Convert TSV to .hic format
    juicer_tools pre \
      -r $(join_by , $(seq 1000 1000 20000)) \
      ${sub_dir}/${sample}-${region}.pre.tsv \
      ${sub_dir}/${sample}-${region}.hic \
      hg38

    # Generate .hic file of summed replicate
    if [[ ${sample} == *'-2' ]]; then
       juicer_tools pre \
         -r $(join_by , $(seq 1000 1000 20000)) \
         <(cat  ${sub_dir}/${sample/-*/}*-${region}.pre.tsv) \
         ${sub_dir}/${sample/-*/}-${region}-sum.hic \
         hg38
    fi

    cp *.hic /media/stephen/Data/hic_01/data/hic_formatted/

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

  rm "${data_dir}"/"${sample}".R[12].filt.bam
  rm "${diffhic_dir}"/"${sample}".captured.sam

done

# Run hicexplorer script
while IFS=$'\t' read -r chr start end region; do
  for binsize in $(seq 1000 1000 10000); do
    /home/stephen/phd/scripts/hicexplorer_normalize_v2.sh \
      -r "${region}" \
      -c "${chr}" \
      -s "${start}" \
      -e "${end}" \
      -b "${binsize}" \
      -t 4 \
      -d "${hcx_dir}"/all_regions/"${region}" "${samples[@]}"
  done
done <"${capture_regions}"

# Switch to py27 environment for hifive
eval "$(conda shell.bash hook)"
conda activate py27env

# Run QUASAR-QC on all regions at multiple bin ranges

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

