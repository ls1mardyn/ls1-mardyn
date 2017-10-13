#!/bin/bash
#@ job_type = parallel
#@ energy_policy_tag = NONE
#@ initialdir=/home/hpc/pr48te/di56giq5/svn/WR_2017/SuperMUC-Phase1/particleTest-8-nodes
#@ job_name = LLRUN
#@ class = test
#@ node_usage = not_shared
#@ wall_clock_limit = 00:30:00
#@ network.MPI = sn_all,not_shared,us
#@ notification = error
#@ notify_user = tchipev@in.tum.de
#@ output = jobout/LLRUN.out.$(jobid)
#@ error =  jobout/LLRUN.err.$(jobid)
#@ node = 8
#@ island_count=1,1
#@ total_tasks = 8
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
export OMP_NUM_THREADS=32
# Pinning, this will use the physical and virtual cores
export MP_TASK_AFFINITY=cpu:$OMP_NUM_THREADS
mpiexec -n 8 ./MarDyn_5598.PAR_RELEASE_AVX 8max.xml --steps 6 --final-checkpoint=0
