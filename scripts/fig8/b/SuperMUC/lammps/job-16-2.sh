#!/bin/bash
#@ energy_policy_tag = tchipev_lammps
#@ minimize_time_to_solution = yes
##
## optional: energy policy tags
#
# DO NOT USE environment = COPY_ALL
#@ job_type = MPICH
#@ class = test
#@ node = 1
#@ island_count=1,1
#@ total_tasks=16
#@ wall_clock_limit = 00:30:00
#@ job_name = lammps-first
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/svn/WR_2017/SuperMUC-Phase1/lammps-comparison/lammps-final-their-lj-newton-on/
#@ output = output-of-jobs/job$(jobid).out
#@ error = output-of-jobs/job$(jobid).err
#@ notification=always
#@ notify_user=tchipev@in.tum.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
#setup of environment
module unload mpi.ibm
module load mpi.intel
module switch intel intel/18.0
module switch mkl mkl/2018
module load tbb/2018

#0000000 from USER-INTEL/TESTS/README
export I_MPI_PIN_DOMAIN=core
export I_MPI_FABRICS=shm		# For single node

export OMP_NUM_THREADS=2

suffs="nocheck check check-big"
rcs="rc25 rc35 rc50"

for i in {3..7} ;
do
	for rc in $rcs ;
	do
		for suf in $suffs ; 
		do
			inputfilename="in.intel.lj-$suf-$rc"
			outputfilename="output-of-jobs/out-lammps-16-2-$inputfilename-$i.txt"
			mpiexec -n 16 ../lmp_intel_cpu_intelmpi -in $inputfilename -log none > $outputfilename
		done
	done
done
