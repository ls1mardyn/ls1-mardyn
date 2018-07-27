#!/bin/bash
#@ job_type = parallel
#@ energy_policy_tag = NONE
##@ group = vip
##@ set_frequency = 2.3
#@ initialdir=/home/hpc/pr48te/di25hup2/2017-10-21-benzene/scaling/real_runs/strongc08_2017-11-9-noglobaloverlapp
#@ job_name = mardyn-strong-128-benzen-c08
#@ class = general
#@ node_usage = not_shared
#@ wall_clock_limit = 1000
#@ network.MPI = sn_all,not_shared,us
#@ notification = error
#@ notify_user = seckler@in.tum.de
#@ output = output-of-jobs/LLRUN.out.128.$(jobid)
#@ error =  output-of-jobs/LLRUN.err.128.$(jobid)
#@ node = 128
#@ island_count = 1,1
#@ total_tasks = 256
##@ hold = user
#@ queue



### @ Other people using this script - change the Email address line! notify_user = xx@in.tum.de

NumProcs=256
#inputFileName=128-nodes.xml
inputFileName=config.xml
executableName=MarDyn-gcc72-ibmmpi-14-DOUBLE-5713-nooverlapping


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

envfilename=output-of-jobs/LLRUN.env.128
testhpfilename=output-of-jobs/LLRUN.thp.128

echo "env var before"
echo " " >>$envfilename
echo " " >>$envfilename
echo " " >>$envfilename
env >> $envfilename

export OMP_NUM_THREADS=16
export MP_TASK_AFFINITY=cpu:$OMP_NUM_THREADS

export MP_SINGLE_THREAD=no
export MP_USE_BULK_XFER=yes
export MP_CSS_INTERRUPT=yes

# should not be needed for strong scaling, ge?
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
