#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 8:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 10   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "mnase R"   # job name
#SBATCH --mem=100G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=job.slurm.txt
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hugj2006@iastate.edu   # email address
#SBATCH --partition=whatever -- use sinfo to get queue names


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

# module load r/3.5.0-py2-ufvuwmm
# R CMD BATCH 5a.WGCNA.r
bash runiSeg.sh