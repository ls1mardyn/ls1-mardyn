#!/bin/bash
#@ energy_policy_tag = tchipev_lammps
#@ minimize_time_to_solution = yes
#@ job_type = MPICH
#@ class = test
#@ node = 1
#@ island_count=1,1
#@ total_tasks=1
#@ wall_clock_limit = 00:30:00
#@ job_name = lammps-first
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/svn/WR_2017/SuperMUC-Phase1/lammps-comparison/water/2-lammps
#@ output = out/job$(jobid).out
#@ error = out/job$(jobid).err
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

export MPI_NUM_PROCESS=1
export OMP_NUM_THREADS=32

#0000000 from USER-INTEL/TESTS/README
export I_MPI_PIN_DOMAIN=omp
export I_MPI_FABRICS=shm		# For single node


for i in {0..3} ;
do
	mpiexec -n $MPI_NUM_PROCESS ../../lmp_intel_cpu_intelmpi -sf omp -pk omp $OMP_NUM_THREADS -in tip4p.in -log out/log-$MPI_NUM_PROCESS-$OMP_NUM_THREADS-$i.txt
done
