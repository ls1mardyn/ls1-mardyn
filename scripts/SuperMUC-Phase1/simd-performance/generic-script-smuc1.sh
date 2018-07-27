#!/bin/bash
#@ job_type = parallel
#@ energy_policy_tag = tchipev00
#@ minimize_time_to_solution = yes
###@ ll_res_id=srv03-ib.56.r
###@ energy_policy_tag = lu78toq.tag_min_t2s
###@ group = vip
###@ set_frequency = 2.3
#@ initialdir=/home/hpc/pr48te/di56giq5/svn/WR_2017/SuperMUC-Phase1/simd-performance
#@ job_name = GENERIC_ARG_NAME
#@ class = GENERIC_ARG_CLASS
#@ node_usage = not_shared
#@ wall_clock_limit = GENERIC_ARG_TIME
#@ network.MPI = sn_all,not_shared,us
#@ notification = always
#@ notify_user = tchipev@in.tum.de
#@ output = output-of-jobs-GENERIC_ARG_NODES-nodes/LLRUN.out.GENERIC_ARG_TOTAL_TASKS.27Ghz.$(jobid)
#@ error =  output-of-jobs-GENERIC_ARG_NODES-nodes/LLRUN.err.GENERIC_ARG_TOTAL_TASKS.27Ghz.$(jobid)
#@ node = GENERIC_ARG_NODES
#@ island_count=GENERIC_ARG_ISLAND_COUNT
#@ total_tasks = GENERIC_ARG_TOTAL_TASKS
#@ queue



### @ Other people using this script - change the Email address line! notify_user = xx@in.tum.de

NumProcs=GENERIC_ARG_TOTAL_TASKS
executableName=GENERIC_ARG_EXECUTABLE


#-------------------------------------
#Settings for POE
#-------------------------------------
. /etc/profile
. /etc/profile.d/modules.sh

module unload mpi.ibm
module load mpi.ibm/1.4_gcc


export MP_SINGLE_THREAD=no
export OMP_NUM_THREADS=GENERIC_ARG_NUMOMP

export MP_TASK_AFFINITY=cpu:$OMP_NUM_THREADS

export MP_SINGLE_THREAD=no
export MP_USE_BULK_XFER=no

iVec=GENERIC_ARG_VEC
iRc=GENERIC_ARG_RC
iRepe=GENERIC_ARG_REPE
iSize=GENERIC_ARG_SIZE
schemes="slice c08"
echo "scheme: $iScheme"
inputFileName=lj-rc$iRc.xml
outputFileName="output-of-jobs-GENERIC_ARG_NODES-nodes/out-$iVec-$iRc-$iRepe.txt"
./$executableName $inputFileName --steps 9 --final-checkpoint=0 >$outputFileName
