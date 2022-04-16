#!/bin/bash
#  File gscript.pbs
# ~Dec. 2021 gtayana@rt.sdstate.edu
# Runs a single node R job on a node of the cluster.

#SBATCH --job-name=BME
#SBATCH --ntasks-per-node=40
#SBATCH --time=7-00:00:00
#SBATCH --partition=bigmem
#SBATCH --output=out-%j-%N.log  # Standard output (log file)

module load R/4.0.2-gcc-6.5.0-mkl 

date
time R --no-save < Multi_trait_one_Env_BME2015.R &> gsnerfresult_out.log
date
