#!/bin/bash
#@ energy_policy_tag = tchipev00
#@ minimize_time_to_solution = yes
#@ job_type = parallel
#@ initialdir = $(home)/svn/WR_2017/SuperMUC-Phase1/lammps-comparison/ls1mardyn-their-lj
#@ class = test
#@ node_usage = not_shared
#@ wall_clock_limit = 00:30:00
#@ job_name = ls1-lam-32-1
#@ network.MPI = sn_all,not_shared,us
#@ notification = always
#@ notify_user = tchipev@in.tum.de
#@ output = output-of-jobs/job$(jobid).out
#@ error = output-of-jobs/job$(jobid).err
#@ node = 1
#@ island_count=1,1
#@ total_tasks = 32
#@ queue



### @ Other people using this script - change the Email address line! notify_user = xx@in.tum.de

executableName=MarDyn_5713.PAR_RELEASE_AVX-gcc-7-ibmmpi-RMM-SINGLE

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

#export MP_SINGLE_THREAD=no
#export MP_USE_BULK_XFER=no

Schemes="slice c08"
for iScheme in $Schemes ;
do
	inputfilename=1-nodes-$iScheme.xml
	for i in {0..3} ;
	do
		outputFileName=output-of-jobs/out-ls1-32-1-$iScheme-$i.txt
		mpiexec -n 32 ../$executableName $inputfilename --steps 101 --final-checkpoint=0 >$outputFileName
	done
done
