#!/bin/bash
#SBATCH -J DC_GV
#SBATCH -o ../../../log/%x.out
#SBATCH -D ../../../build/src
#SBATCH --partition=medium
#SBATCH --get-user-env
#SBATCH --mail-type=all
#SBATCH --mem=800mb
#SBATCH --mail-user=hocksa@hsu-hh.de
#SBATCH --export=NONE
#SBATCH --time=00:10:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18
#SBATCH --error=../../../log/%x.err

export OMP_NUM_THREADS=18
ulimit -s 1000000

module load slurm_setup
module load mpi/2021.6.0

mpiexec -n $((node*4)) ./MarDyn ../../examples/AdResS/DropletCoalescence/vap/config_3_generateVap.xml --loop-abort-time=570