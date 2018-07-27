#!/bin/bash
#@ job_type = parallel
#@ ll_res_id=srv03-ib.56.r
#@ energy_policy_tag = lu78toq.tag_min_t2s
#@ group = vip
#@ set_frequency = 2.7
#@ initialdir=/home/hpc/pr48te/di56giq5/svn/WR_2017/SuperMUC-Phase1/2017-10-24/1clj-double-rc35/weak
#@ job_name = GENERIC_ARG_NAME
#@ class = GENERIC_ARG_CLASS
#@ node_usage = not_shared
#@ wall_clock_limit = GENERIC_ARG_TIME
#@ network.MPI = sn_all,not_shared,us
#@ notification = always
#@ notify_user = tchipev@in.tum.de
#@ output = output-of-jobs/LLRUN.out.GENERIC_ARG_TOTAL_TASKS.$(jobid)
#@ error =  output-of-jobs/LLRUN.err.GENERIC_ARG_TOTAL_TASKS.$(jobid)
#@ node = GENERIC_ARG_NODES
#@ island_count=GENERIC_ARG_ISLAND_COUNT
#@ total_tasks = GENERIC_ARG_TOTAL_TASKS
#@ hold=user
#@ queue



### @ Other people using this script - change the Email address line! notify_user = xx@in.tum.de

NumProcs=GENERIC_ARG_TOTAL_TASKS
inputFileName=GENERIC_ARG_NODES-nodes.xml
executableName=GENERIC_ARG_EXECUTABLE


#-------------------------------------
#Settings for POE
#-------------------------------------
. /etc/profile
. /etc/profile.d/modules.sh

module unload mpi.ibm
module load mpi.ibm/1.4_gcc

envfilename=output-of-jobs/LLRUN.env.GENERIC_ARG_TOTAL_TASKS
testhpfilename=output-of-jobs/LLRUN.thp.GENERIC_ARG_TOTAL_TASKS

echo "env var before"
echo " " >>$envfilename
echo " " >>$envfilename
echo " " >>$envfilename
env >> $envfilename

export MP_SINGLE_THREAD=no
export OMP_NUM_THREADS=32

export MP_TASK_AFFINITY=cpu:$OMP_NUM_THREADS

export MP_SINGLE_THREAD=no
export MP_USE_BULK_XFER=no
export MP_BUFFER_MEM=8M,32M
export MP_EAGER_LIMIT=64
export MP_RFIFO_SIZE=1048576


echo "env var after"
echo " " >>$envfilename
echo " " >>$envfilename
echo " " >>$envfilename
env >> $envfilename


../../testhp.sh > $testhpfilename &

mpiexec -n $NumProcs ../../$executableName $inputFileName --steps 11 --final-checkpoint=0
