#!/bin/bash

#SBATCH -o /home/hpc/pr63so/lu32reb2/trunk/examples/scenarios/homogeneous-1clj/legacy/16.out
#SBATCH -D /home/hpc/pr63so/lu32reb2/trunk/examples/scenarios/homogeneous-1clj/legacy/
#SBATCH -J h1clj-016
#SBATCH --get-user-env
#SBATCH --partition=snb
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=all
#SBATCH --mail-user=eckhardw@in.tum.de
#SBATCH --export=NONE
#SBATCH --time=01:00:00

source /etc/profile.d/modules.sh

export NTIMES=100

# -ppn: processes per node  -n: #processes
mpiexec.hydra -ppn 16 -n 16 ../../../../src/MarDyn.PAR_RELEASE ljfluid_640k_rc3.cfg $NTIMES --final-checkpoint=0 
