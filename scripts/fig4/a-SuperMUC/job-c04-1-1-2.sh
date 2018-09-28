#!/bin/bash
#@ job_type = parallel
#@ energy_policy_tag = tchipev00
#@ minimize_time_to_solution = yes
###@ ll_res_id=srv03-ib.56.r
###@ energy_policy_tag = lu78toq.tag_min_t2s
###@ group = vip
###@ set_frequency = 2.3
#@ initialdir=/home/hpc/pr48te/di56giq5/svn/WR_2017/SuperMUC-Phase1/omp-strong-2/
#@ job_name = omp-strong-1-2-sp-35
#@ class = test
#@ node_usage = not_shared
#@ wall_clock_limit = 00:30:00
#@ network.MPI = sn_all,not_shared,us
#@ notification = always
#@ notify_user = tchipev@in.tum.de
#@ output = output-of-jobs-1-nodes/LLRUN.out.1.27Ghz.$(jobid)
#@ error =  output-of-jobs-1-nodes/LLRUN.err.1.27Ghz.$(jobid)
#@ node = 1
#@ island_count=1,1
#@ total_tasks = 1
#@ queue



### @ Other people using this script - change the Email address line! notify_user = xx@in.tum.de

NumProcs=1
executableName=MarDyn-5713-patched-with-c04-RMM-single-OpenMP-noMPI


#-------------------------------------
#Settings for POE
#-------------------------------------
. /etc/profile
. /etc/profile.d/modules.sh

module unload mpi.ibm
module load mpi.ibm/1.4_gcc


export MP_SINGLE_THREAD=no
export OMP_NUM_THREADS=2

export MP_TASK_AFFINITY=cpu:$OMP_NUM_THREADS

export MP_SINGLE_THREAD=no
export MP_USE_BULK_XFER=no

# for loop over system sizes
# ONE repetition for now
#iRepe=GENERIC_ARG_REPE
for ((iRepe=0; iRepe <= 4; iRepe +=1)) ; 
do
	echo "repetition $iRepe"
	iSize=32
	#for ((iSize=4096; iSize >= 1; iSize /= 2)) ;
	#do
		echo "system size: $iSize"
		#for loop over schemes
		schemes="c04"
		for iScheme in $schemes ;
		do
			echo "scheme: $iScheme"
			inputFileName=1-nodes-$iSize-$iScheme.xml
			outputFileName="output-of-jobs-1-nodes/out-$NumProcs-$OMP_NUM_THREADS-$iSize-$iScheme-$iRepe.txt"
			mpiexec -n $NumProcs ./$executableName $inputFileName --steps 11 --final-checkpoint=0 >$outputFileName
		done
	#done
done
