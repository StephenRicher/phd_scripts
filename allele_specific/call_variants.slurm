#!/usr/bin/env bash

#SBATCH --account=ba-fs6sdb
#SBATCH --array=1-3
#SBATCH --job-name=call_variants_%A_%a
#SBATCH --output=/home/u/sr467/scratch/projects/hic-01/qc/call_variants_%A_%a.out
#SBATCH --error=/home/u/sr467/scratch/projects/hic-1/qc/call_variants_%A_%a.err
#SBATCH --partition=batch-128gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr467@bath.ac.uk

module purge
module load slurm
module load openmpi/intel
module load java/jdk/1.8.0
module load bcftools/1.9
module load samtools/1.9

# Excluding replicate number
samples=( 'HB2_CL4' 'HB2_WT' 'MCF7' )

# File containing paths to processed BAMS for variant calling.
processed_bams="/home/u/sr467/scratch/projects/hic-01/paths/hic-processed_bams.txt"

# Get sample-specific BAMs from processed BAMS file
index=$(("${SLURM_ARRAY_TASK_ID}"-1))
sample="${samples["${index}"]}"
input_bams=( $(grep "${sample}" "${process_bams}") )

# Reference genome (unmasked) that processed BAMs were mapped to.
genome="/home/u/sr467/scratch/projects/genomes/GRCh38/wgs/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Capture regions from capture HiC protocol
capture_regions="/home/u/sr467/scratch/scripts/phd_scripts/capture_regions.bed"

# Other options
data_dir="/home/u/sr467/scratch/projects/hic-01/allele_specific"
threads="${SLURM_CPUS_PER_TASK}"

# Run call_variants script
mpirun call_variants \
    -g "${genome}" \
    -c "${capture_regions}" \
    -d "${data_dir}" \
    -j "${threads}" \
    "${input_bams[@]}"
