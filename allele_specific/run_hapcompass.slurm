#!/usr/bin/env bash

#SBATCH --account=free
#SBATCH --array=1-63
#SBATCH --job-name=run_hapcompass_%A_%a
#SBATCH --output=/home/u/sr467/scratch/projects/hic-01/qc/run_hapcompass_%A_%a.out
#SBATCH --error=/home/u/sr467/scratch/projects/hic-01/qc/run_hapcompass_%A_%a.err
#SBATCH --partition=batch-128gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr467@bath.ac.uk

module purge
module load slurm
module load openmpi/intel
module load java/jdk/1.8.0
module load bcftools/1.9

# Excluding replicate number
samples=( 'HB2_CL4' 'HB2_WT' 'MCF7' )

# Capture regions from capture HiC protocol
capture_regions="/home/u/sr467/scratch/scripts/phd_scripts/capture_regions.bed"

# Other options
data_dir="/home/u/sr467/scratch/projects/hic-01/allele_specific"
threads="${SLURM_CPUS_PER_TASK}"

count=1
while IFS=$'\t' read -r chr start end region; do
    for sample in "${samples[@]}"; do
        if [[ "${count}" == "${SLURM_ARRAY_TASK_ID}" ]]; then
            mpirun ./run_hapcompass.sh \
                -v "${data_dir}"/"${sample}".sorted.vcf.gz \
                -c "${chr}" \
                -s "${start}" \
                -e "${end}" \
                -r "${region}" \
                -d "${data_dir}" \
                -j "${threads}"
        fi
        ((count=count+1))
    done
done <"${capture_regions}"
