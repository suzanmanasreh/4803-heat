#!/usr/bin/env bash

# change my email to yours
# change the output file path to something different so you don't overwrite it

#SBATCH -N3 --ntasks-per-node=24
#SBATCH -t1:00:00
#SBATCH -q coc-ice
#SBATCH --mail-user=smanasreh6@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "./run_logs/run.log"

make clean
make par

mkdir 2dp

ml openmpi boost
srun -n 64 ./main_2d_par.o