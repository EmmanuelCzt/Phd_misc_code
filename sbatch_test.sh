#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=test_sbatch

### Limit run time "days-hours:minutes:seconds"
#SBATCH --time=01:00:00

### Requirements
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8GB
#SBATCH --account=primate_hic

### Email
#SBATCH --mail-user=emmanuel.cazottes@etu.u-paris.fr
#SBATCH --mail-type=ALL

### Output
#SBATCH --output=/shared/home/<user>/demojob-%j.out

################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
echo '########################################'

# modules loading
module samtools


# What you actually want to launch
ASSEMBLY=$("hg38")
INPATH=$(primates_multiple_alignment/assemblies/)
OUTPATH=$(primates_multiple_alignment/primaX/)

mkdir ${OUTPATH}

samtools faix ${INPATH}${ASSEMBLY}.fa chrX > ${OUTPATH}${ASSEMBLY}.chrX.fa

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
