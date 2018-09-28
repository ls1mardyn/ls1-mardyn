#!/bin/bash
#@ job_type = parallel
#@ energy_policy_tag = tchipev00
#@ minimize_time_to_solution = yes
###@ ll_res_id=srv03-ib.56.r
###@ energy_policy_tag = lu78toq.tag_min_t2s
###@ group = vip
###@ set_frequency = 2.3
#@ initialdir=/home/hpc/pr48te/di56giq5/svn/WR_2017/SuperMUC-Phase1/omp-strong-2/
#@ job_name = mpi-strong-16-32-sp-35
#@ class = test
#@ node_usage = not_shared
#@ wall_clock_limit = 00:30:00
#@ network.MPI = sn_all,not_shared,us
#@ notification = always
#@ notify_user = tchipev@in.tum.de
#@ output = output-of-jobs-mpi/LLRUN.out.16.27Ghz.$(jobid)
#@ error =  output-of-jobs-mpi/LLRUN.err.16.27Ghz.$(jobid)
#@ node = 1
#@ island_count=1,1
#@ total_tasks = 16
#@ queue



### @ Other people using this script - change the Email address line! notify_user = xx@in.tum.de

NumProcs=16
executableName=MarDyn_5713.PAR_RELEASE_AVX-gcc-7-ibmmpi-RMM-SINGLE-no-openmp


#-------------------------------------
#Settings for POE
#-------------------------------------
. /etc/profile
. /etc/profile.d/modules.sh

module unload mpi.ibm
module load mpi.ibm/1.4_gcc


export MP_SINGLE_THREAD=no
export OMP_NUM_THREADS=1

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
		#schemes="slice c08"
		schemes="slice"
		for iScheme in $schemes ;
		do
			echo "scheme: $iScheme"
			inputFileName=1-nodes-$iSize-$iScheme.xml
			outputFileName="output-of-jobs-mpi/out-$NumProcs-$OMP_NUM_THREADS-$iSize-$iScheme-$iRepe.txt"
			mpiexec -n $NumProcs ./$executableName $inputFileName --steps 11 --final-checkpoint=0 >$outputFileName
		done
	#done
done
