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
	"hictools truncate --restriction ^GATC -zu \
    	--summary "${qc}"/hictools-truncate_summary.txt {1} \
    	> {=1 s/trim/trim-trunc/ =}" \
	::: "${data_dir}"/*trim.fq.gz

## Run FastQC on adapter trimmed and truncated data ##
fastqc --threads 12 --outdir "${qc}" "${data_dir}"/*trim-trunc.fq.gz

# Map R1 and R2 reads
for sample in "${samples[@]}"; do
	hictools map --index "${genome_index}" --sample "${sample}" \
		--bowtie2_log "${qc}"/"${sample}"_alignment_stats.txt \
		--log "${qc}"/"${sample}".bowtie2.logfile \
		--output "${data_dir}"/"${sample}".bam \
		"${data_dir}"/"${sample}"-R[14]-trim-trunc.fq.gz

	hictools deduplicate --sample "${sample}" \
		--output "${data_dir}"/"${sample}".dedup.bam \
		--log "${qc}"/"${sample}".dedup.logfile \
		"${data_dir}"/"${sample}".bam

	rm "${data_dir}"/"${sample}".bam

	hictools process --digest "${genome_digest}" \
		--log "${qc}"/"${sample}".process.logfile \
		--output "${data_dir}"/"${sample}".proc.bam \
		"${data_dir}"/"${sample}".dedup.bam
	
	rm "${data_dir}"/"${sample}".dedup.bam

	hictools extract --sample "${sample}" --gzip \
		--log "${qc}"/"${sample}".extract.logfile
		<(samtools view -s 42.1 "${data_dir}"/"${sample}".proc.bam) \
		>> "${qc}"/hic_filter_qc.txt

	# CONFIRM THESE PARAMETERS!!!##
	hictools filter --max_ditag 1000 --min_outward 1000 \
		--log "${qc}"/"${sample}".filter.logfile \
		--output "${data_dir}"/"${sample}".filt.bam \
		"${data_dir}"/"${sample}".proc.bam

	rm "${data_dir}"/"${sample}".proc.bam
done
