#!/bin/bash
#@ job_type = parallel
#@ energy_policy_tag = NONE
#@ initialdir=/home/hpc/pr48te/di56giq5/svn/WR_2017/SuperMUC-Phase1/benzene/sim01/run02
#@ job_name = LLRUN
#@ class = test
#@ node_usage = not_shared
#@ wall_clock_limit = 00:30:00
#@ network.MPI = sn_all,not_shared,us
#@ notification = error
#@ notify_user = tchipev@in.tum.de
#@ output = output-of-jobs/LLRUN.out.$(jobid)
#@ error =  output-of-jobs/LLRUN.err.$(jobid)
#@ node = 1
#@ island_count=1,1
#@ total_tasks = 1
#@ queue



### @ Other people using this script - change the Email address line! notify_user = xx@in.tum.de



#-------------------------------------
#Settings for POE
#-------------------------------------
. /etc/profile
. /etc/profile.d/modules.sh
module unload mpi.ibm
module load mpi.ibm/1.4_gcc
export MP_SINGLE_THREAD=no
export OMP_NUM_THREADS=1
# Pinning, this will use the physical and virtual cores
export MP_TASK_AFFINITY=cpu:$OMP_NUM_THREADS
mpiexec -n 1 ./MarDyn_5670_-gcc-7-same-NONRMM config.xml --steps 11 --final-checkpoint=0
