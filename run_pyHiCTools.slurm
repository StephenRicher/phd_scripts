#!/bin/bash

#SBATCH --account=ba-fs6sdb
#SBATCH --array=1-6
#SBATCH --job-name=run_pyHiCTools_%A_%a
#SBATCH --output=/home/u/sr467/scratch/projects/hic-01/job_logs/run_pyHiCTools_%A_%a.out
#SBATCH --error=/home/u/sr467/scratch/projects/hic-01/job_logs/run_pyHiCTools_%A_%a.err
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --workdir=/home/u/sr467/scratch/projects/hic-01/
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr467@bath.ac.uk

module purge
module load slurm
module load openmpi/intel
module load samtools/1.9
module load bowtie2
module load python3

# File containing paths to trimmed FASTQ files for mapping
trimmed_fastqs="/home/u/sr467/scratch/projects/hic-01/paths/trimmed_fastqs.txt"

# Bowtie2 index
bt2_idx=/home/u/sr467/scratch/projects/genomes/GRCh38/bt2_index/GRCh38

# Digest file from pyHiCTools digest
digest=/home/u/sr467/scratch/projects/genomes/GRCh38/GRCh38_Mbo1-digest.txt.gz

# Restriction enzyme cut sequence
re_seq='^GATC'

# Other options
data_dir="/home/u/sr467/scratch/projects/hic-01/data"
qc_dir="/home/u/sr467/scratch/projects/hic-01/qc"
threads="${SLURM_CPUS_PER_TASK}"

# Get nth and nth-1 line from trimmed_fastqs
line_num=$( bc <<< "2 * "${SLURM_ARRAY_TASK_ID}"" )
input_reads=( $(head -n "${line_num}" "${trimmed_fastqs}" | tail -n 2) )

# Run pyHiCTools
mpirun run_pyHiCTools \
    -f "${input_reads[0]}" \
    -r "${input_reads[1]}" \
    -x "${bt2_idx}" \
    -i "${digest}" \
    -s "${re_seq}" \
    -d "${data_dir}" \
    -q "${qc_dir}" \
    -j "${threads}"


