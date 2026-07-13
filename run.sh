#!/bin/bash
#SBATCH --job-name="lammmps HTC"
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G

module load 2024r1
module load openmpi/4.1.6
module load lammps/20230802

lmp -in input.in > input.out
