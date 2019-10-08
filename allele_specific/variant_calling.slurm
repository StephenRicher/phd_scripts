#!/usr/bin/env bash

#SBATCH --account=ba-fs6sdb
#SBATCH --job-name=variant_calling_%j
#SBATCH --output=/home/u/sr467/scratch/projects/HiC/allele_specific/variant_calling_%j.out
#SBATCH --error=/home/u/sr467/scratch/projects/HiC/allele_specific/variant_calling_%j.err
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

sample=${1}
shift

outdir="/home/u/sr467/scratch/projects/HiC/allele_specific"
genome="/home/u/sr467/scratch/projects/genomes/GRCh38/wgs/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
regions="/home/u/sr467/scratch/scripts/phd_scripts/capture_regions.bed"
hapcompass="/home/u/sr467/scratch/scripts/phd_scripts/allele_specific/hapcompass_v0.8.2/hapcompass.jar"
hc2vcf="/home/u/sr467/scratch/scripts/phd_scripts/allele_specific/hapcompass_v0.8.2/hc2vcf.jar"

mpirun /home/u/sr467/scratch/scripts/phd_scripts/allele_specific/variant_calling.sh -s ${sample} -o ${outdir} -g ${genome} -r ${regions} -h ${hapcompass} -v ${hc2vcf} -t ${SLURM_CPUS_PER_TASK} "${@}"
