#!/bin/bash
#@ job_type = parallel
##@ energy_policy_tag = lu78toq.tag_min_t2s
#@ energy_policy_tag = NONE
##@ group = vip
##@ set_frequency = 2.3
#@ ll_res_id = GENERIC_ARG_RESERVATION
#@ initialdir=/home/hpc/pr48te/di25hup2/2017-10-21-benzene/scaling/real_runs/weakc08-bigger
#@ job_name = GENERIC_ARG_NAME
#@ class = GENERIC_ARG_CLASS
#@ node_usage = not_shared
#@ wall_clock_limit = GENERIC_ARG_TIME
#@ network.MPI = sn_all,not_shared,us
#@ notification = error
#@ notify_user = seckler@in.tum.de
#@ output = output-of-jobs/LLRUN.out.GENERIC_ARG_NODES.$(jobid)
#@ error =  output-of-jobs/LLRUN.err.GENERIC_ARG_NODES.$(jobid)
#@ node = GENERIC_ARG_NODES
#@ island_count = GENERIC_ARG_ISLAND_COUNT
#@ total_tasks = GENERIC_ARG_TOTAL_TASKS
##@ hold = user
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
#module unload mpi.intel
module load mpi.ibm/1.4_gcc
#module load mpi.intel/2017
module switch gcc gcc/7

envfilename=output-of-jobs/LLRUN.env.GENERIC_ARG_NODES
testhpfilename=output-of-jobs/LLRUN.thp.GENERIC_ARG_NODES

echo "env var before"
echo " " >>$envfilename
echo " " >>$envfilename
echo " " >>$envfilename
env >> $envfilename

export MP_SINGLE_THREAD=no
export OMP_NUM_THREADS=16

export MP_TASK_AFFINITY=cpu:$OMP_NUM_THREADS

export MP_SINGLE_THREAD=no
export MP_USE_BULK_XFER=no
# not needed, as we do not really fill the nodes!
#export MP_BUFFER_MEM=8M,32M
#export MP_EAGER_LIMIT=64
#export MP_RFIFO_SIZE=1048576


echo "env var after"
echo " " >>$envfilename
echo " " >>$envfilename
echo " " >>$envfilename
env >> $envfilename


./../testhp.sh > $testhpfilename &

mpiexec -n $NumProcs /home/hpc/pr48te/di25hup2/2017-10-21-benzene/executables/$executableName $inputFileName --steps 11 --final-checkpoint=0
