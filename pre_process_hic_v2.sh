#!/usr/bin/env bash

shopt -s extglob
export TMPDIR=/home/stephen/x_db/DBuck/s_richer/tmp/

## Function to record system memory usage ## 
memlog() {
	echo "      date     time $(free -m | grep total | sed -E 's/^    (.*)/\1/g')"
	while true; do
    echo "$(date '+%Y-%m-%d %H:%M:%S') $(free -m | grep Mem: | sed 's/Mem://g')"
    sleep 1
	done
}

## Set project directory ##
project=/home/stephen/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis2
qc="${project}"/qc
mkdir -p "${qc}"

# TEMPORARILY MODIF PROJECT - THIS WILL BE UPDATED WHEN FILES ARE TRANSFERRED"
proj=/home/stephen/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis
data_dir="${proj}"/raw_data

samples=( HB2_CL4-1 HB2_CL4-2 HB2_WT-1 HB2_WT-2 MCF7-1 MCF7-2 )

## Begin recording memory usage ##
memlog > "${qc}"/memory_usage.logfile &

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
fastq_screen --aligner bowtie2 --threads 6 --outdir "${qc}" "${raw_data[@]}"

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

alignment_stats=alignment_stats.txt
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

	total_pairs=$(samtools view -cf 64 "${intermediate}")
	both_pair_unmapped=$(samtools view -cf 12 "${intermediate}")
	r1_map_r2_unmap=$(samtools view -c -f 72 -F 4 "${intermediate}")
	r2_map_r1_unmap=$(samtools view -c -f 136 -F 4 "${intermediate}")
	unmapped_pairs=$(( both_pair_unmapped + r1_map_r2_ummap + r2_map_r1_unmap ))
	duplicate_pairs=$(sed -n '4p' "${qc}"/"${sample}".dedup_stats.txt | awk '{print $3/2}')
	total_kept_pairs=$(samtools view -cf 64 "${data_dir}"/"${sample}".proc.bam)
	qual_filtered_pairs=$(( total_pairs - total_kept_pairs - duplicate_pairs - unmapped_pairs ))

	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
		"${sample}" "${total_pairs}" "${both_pair_unmapped}" \
		"${r1_map_r2_unmap}" "${r2_map_r1_unmap}" \
		"${duplicate_pairs}" "${qual_filtered_pairs}" \
		"${total_kept_pairs}" >> "${alignment_stats}"

	hictools extract \
		--sample "${sample}" --gzip \
		--log "${qc}"/"${sample}".extract.logfile \
		<(samtools view -s 42.1 "${data_dir}"/"${sample}".proc.bam) \
		>> "${qc}"/hic_filter_qc.txt

	# CONFIRM THESE PARAMETERS!!!##
	hictools filter \
		--max_ditag 1000 --min_outward 1000 \
		--log "${qc}"/"${sample}".filter.logfile \
		"${data_dir}"/"${sample}".proc.bam
		> "${data_dir}"/"${sample}".filt.bam


done
