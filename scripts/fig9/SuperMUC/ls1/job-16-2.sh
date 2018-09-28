#!/bin/bash
#@ energy_policy_tag = tchipev00
#@ minimize_time_to_solution = yes
#@ job_type = parallel
#@ initialdir = $(home)/svn/WR_2017/SuperMUC-Phase1/lammps-comparison/water/mardyn
#@ class = test
#@ node_usage = not_shared
#@ wall_clock_limit = 00:30:00
#@ job_name = ls1-lam-16-2
#@ network.MPI = sn_all,not_shared,us
#@ notification = always
#@ notify_user = tchipev@in.tum.de
#@ output = out/job$(jobid).out
#@ error = out/job$(jobid).err
#@ node = 1
#@ island_count=1,1
#@ total_tasks = 16
#@ queue



### @ Other people using this script - change the Email address line! notify_user = xx@in.tum.de

executableName=MarDyn_6153.PAR_RELEASE_AVX

#-------------------------------------
#Settings for POE
#-------------------------------------
. /etc/profile
. /etc/profile.d/modules.sh

module unload mpi.ibm
module load mpi.ibm/1.4_gcc


export MP_SINGLE_THREAD=no
export OMP_NUM_THREADS=2
export MPI_NUM_PROCESS=16

export MP_TASK_AFFINITY=cpu:$OMP_NUM_THREADS

#export MP_SINGLE_THREAD=no
#export MP_USE_BULK_XFER=no

Schemes="sli c08"
for iScheme in $Schemes ;
do
	inputfilename=config-$iScheme.xml
	for i in {0..3} ;
	do
		outputFileName=out/out-ls1-$MPI_NUM_PROCESS-$OMP_NUM_THREADS-$iScheme-$i.txt
		mpiexec -n $MPI_NUM_PROCESS ../../$executableName $inputfilename --steps 101 --final-checkpoint=0 >$outputFileName
	done
done
