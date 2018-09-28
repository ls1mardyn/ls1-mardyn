#!/bin/bash
#@ job_type = parallel
#@ energy_policy_tag = tchipev00
#@ minimize_time_to_solution = yes
#@ initialdir=/home/hpc/pr48te/di56giq5/svn/WR_2017/SuperMUC-Phase1/particleTest-final-8-nodes/rc25
#@ job_name = pTest-8-8-sp-25
#@ class = test
#@ node_usage = not_shared
#@ wall_clock_limit = 00:30:00
#@ network.MPI = sn_all,not_shared,us
#@ notification = always
#@ notify_user = tchipev@in.tum.de
#@ output = output-of-jobs-8-nodes/LLRUN.out.64.27Ghz.$(jobid)
#@ error =  output-of-jobs-8-nodes/LLRUN.err.64.27Ghz.$(jobid)
#@ node = 8
#@ island_count=1,1
#@ total_tasks = 64
#@ queue



### @ Other people using this script - change the Email address line! notify_user = xx@in.tum.de

NumProcs=64
executableName=MarDyn_5713.PAR_RELEASE_AVX-gcc-7-ibmmpi-RMM-SINGLE


#-------------------------------------
#Settings for POE
#-------------------------------------
. /etc/profile
. /etc/profile.d/modules.sh

module unload mpi.ibm
module load mpi.ibm/1.4_gcc


export MP_SINGLE_THREAD=no
export OMP_NUM_THREADS=4

export MP_TASK_AFFINITY=cpu:$OMP_NUM_THREADS

export MP_SINGLE_THREAD=no
export MP_USE_BULK_XFER=no

iRepe=1
echo "repetition: $iRepe"
iScheme=slice
echo "scheme: $iScheme"

# for loop over system sizes
for ((iSize=1024; iSize >= 8; iSize /= 2)) ;
do
	echo "system size: $iSize"
	inputFileName=8-nodes-$iSize-$iScheme.xml
	outputFileName="output-of-jobs-8-nodes/out-$NumProcs-$iSize-$iScheme-$iRepe.txt"
	mpiexec -n $NumProcs ../$executableName $inputFileName --steps 11 --final-checkpoint=0 >$outputFileName
done
