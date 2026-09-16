#!/bin/bash
#PBS -P vk72
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l walltime=08:00:00
#PBS -l wd
#PBS -l storage=scratch/vk72+gdata/vk72

module load nci-parallel/1.0.0a
module load R/4.3.1

export ncores_per_task=1
export ncores_per_numanode=12

mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file cmds.txt --timeout 3600 --status status.txt --output-dir logs

