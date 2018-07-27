#!/bin/bash
#@ job_type = parallel
##@ energy_policy_tag = lu78toq.tag_min_t2s
#@ energy_policy_tag = NONE
##@ group = vip
##@ set_frequency = 2.3
#@ initialdir=/home/hpc/pr48te/di25hup2/2017-10-21-benzene/scaling/real_runs/weakc08-bigger-1lq6q
#@ job_name = mardyn-weak-8-bigger-1lj6q
#@ class = test
#@ node_usage = not_shared
#@ wall_clock_limit = 00:45:00
#@ network.MPI = sn_all,not_shared,us
#@ notification = error
#@ notify_user = seckler@in.tum.de
#@ output = output-of-jobs/LLRUN.out.8.$(jobid)
#@ error =  output-of-jobs/LLRUN.err.8.$(jobid)
#@ node = 8
#@ island_count = 1,1
#@ total_tasks = 16
##@ hold = user
#@ queue



### @ Other people using this script - change the Email address line! notify_user = xx@in.tum.de

NumProcs=16
inputFileName=8-nodes.xml
executableName=MarDyn-gcc72-ibmmpi-14-DOUBLE-5676


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

envfilename=output-of-jobs/LLRUN.env.8
testhpfilename=output-of-jobs/LLRUN.thp.8

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
