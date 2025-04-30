#!/usr/bin/env bash

# change my email to yours
# change the output file path to something different so you don't overwrite it

#SBATCH -N4 --ntasks-per-node=8
#SBATCH -t1:00:00
#SBATCH -q coc-ice
#SBATCH --mail-user=smanasreh6@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "./run_logs/run_omp.log"

ml openmpi boost
make clean
make par

export OMP_NUM_THREADS=8

# mkdir 2dp

# mpirun -np 4 ./main_2d_par.o
srun -n 4 ./main_2d_par.o