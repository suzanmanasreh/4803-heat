#!/usr/bin/env bash

# change my email to yours
# change the output file path to something different so you don't overwrite it

#SBATCH -N36 --ntasks-per-node=8
#SBATCH -t0:05:00
#SBATCH -q coc-ice
#SBATCH --mail-user=smanasreh6@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "./run_logs/run2_36.log"

ml openmpi boost
make clean
make par

export OMP_NUM_THREADS=8

# mkdir 2dp

# mpirun -np 4 ./main_2d_par.o
#!/usr/bin/env bash

make clean
make
cd build/handout

# mpirun -np 1 --cpus-per-proc 8 --bind-to core ./main_2d_par.o

srun -n 36 ./main_2d_par.o
# srun -n 4 ./main_2d_par.o
# srun -n 16 ./main_2d_par.o
# srun -n 4 ./main_2d_par.o
# srun -n 4 ./main_2d_par.o

export OMP_NUM_THREADS=4
echo "OMP_NUM_THREADS=4"
srun -n 36 ./main_2d_par.o

export OMP_NUM_THREADS=2
echo "OMP_NUM_THREADS=2"
srun -n 36 ./main_2d_par.o

export OMP_NUM_THREADS=1
echo "OMP_NUM_THREADS=1"
srun -n 36 ./main_2d_par.o