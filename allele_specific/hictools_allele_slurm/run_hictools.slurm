#!/usr/bin/env bash

#SBATCH --account=ba-fs6sdb
#SBATCH --job-name=variant_calling_%j
#SBATCH --output=/home/u/sr467/scratch/projects/HiC/allele_specific/hic_dir/hictools_allele_%j.out
#SBATCH --error=/home/u/sr467/scratch/projects/HiC/allele_specific/hic_dir/hictools_allele_%j.err
#SBATCH --partition=batch-64gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sr467@bath.ac.uk

module purge
module load slurm
module load openmpi/intel
module load samtools/1.9
module load python3
module load bowtie2

sample=${1}

mpirun /home/u/sr467/scratch/projects/HiC/allele_specific/hic_dir/run_hictools.sh "${sample}"

