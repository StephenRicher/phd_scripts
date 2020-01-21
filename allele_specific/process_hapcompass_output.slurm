#!/usr/bin/env bash

#SBATCH --account=free
#SBATCH --array=1-3
#SBATCH --job-name=call_variants_%A_%a
#SBATCH --output=/home/u/sr467/scratch/projects/hic-01/qc/process_hapcompass_output_%A_%a.out
#SBATCH --error=/home/u/sr467/scratch/projects/hic-1/qc/process_hapcompass_output_%A_%a.err
#SBATCH --partition=batch-all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr467@bath.ac.uk

module purge
module load slurm
module load openmpi/intel
module load java/jdk/1.8.0

### ARCHIVED VERSION DELETE THIS ####

# Excluding replicate number
samples=( 'HB2_CL4' 'HB2_WT' 'MCF7' )

# Get sample-specific BAMs from processed BAMS file
index=$(("${SLURM_ARRAY_TASK_ID}"-1))
sample="${samples["${index}"]}"

# Other options
data_dir="/home/u/sr467/scratch/projects/hic-01/allele_specific"

mpirun run_hapcompass \
    -v "${data_dir}"/"${sample}".sorted.vcf.gz \
    -d "${data_dir}" \
