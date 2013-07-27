#!/bin/bash

#SBATCH -o /home/hpc/pr63so/lu32reb2/MarDyn-jenkins/examples/scenarios/homogeneous-1clj/legacy/8.out
#SBATCH -D /home/hpc/pr63so/lu32reb2/MarDyn-jenkins/examples/scenarios/homogeneous-1clj/legacy/
#SBATCH -J h1clj-008
#SBATCH --get-user-env
#SBATCH --partition=snb
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=all
#SBATCH --mail-user=eckhardw@in.tum.de
#SBATCH --export=NONE
#SBATCH --time=00:09:00

source /etc/profile.d/modules.sh

export NTIMES=1000

# -ppn: processes per node  -n: #processes
mpiexec.hydra -ppn 8 -n 8 ../../../../src/MarDyn.PAR_RELEASE ljfluid_640k_rc3.cfg $NTIMES --final-checkpoint=0 
