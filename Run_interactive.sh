#!/bin/bash -l

#SBATCH --ntasks 28
#SBATCH -J SimulLenS
#SBATCH -o SimulLenS_%a.out
#SBATCH -e SimulLenS_%a.err
#SBATCH -p cosma6
#SBATCH -A dp004
#SBATCH --exclusive
#SBATCH -t 500

module purge
module load intel_comp/2018
module load intel_mpi/2018
module load parallel_hdf5/1.10.3
module load gsl
module load fftw
module load python/2.7.15

model=$1
#model_str=$model'_a'; seed=2080
model_str=$model'_b'; seed=4257

mkdir -p ~/SimulLens/FORGE_lens/node$model_str/kappa/; mkdir -p ~/SimulLens/FORGE_lens/node$model_str/delta/ ; 

mpirun -np 5 ~/array_job_test/parallel_tasks  1 25 "mkdir thread_%d; cp -v che* thread_%d; cp -v CosmoDist.py thread_%d; cp -v SimulLens thread_%d ; cd thread_%d ; ./SimulLens '' %d %d ../../../FORGE/node$model_str/2d_proj/ ../../../FORGE_lens/node$model_str/ ../random_shifts/ $seed  $model; cd .."
#mpirun -np 5 ~/array_job_test/parallel_tasks  1 25 "mkdir thread_%d; cp -v che* thread_%d; cp -v CosmoDist.py thread_%d; cp -v SimulLens_gamma thread_%d ; cd thread_%d ; ./SimulLens_gamma '' %d %d ../../../FORGE/node$model_str/2d_proj/ ../../../FORGE_lens/node$model_str/ ../random_shifts/ $seed  $model; cd .."


