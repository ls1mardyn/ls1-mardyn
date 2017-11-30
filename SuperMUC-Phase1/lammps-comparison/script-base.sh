#!/bin/bash
##
## optional: energy policy tags
#
# DO NOT USE environment = COPY_ALL
#@ job_type = MPICH
#@ class = large
#@ node = 1000
#@ island_count=2,3
#@ total_tasks=8000
## other example
##@ tasks_per_node = 8
#@ wall_clock_limit = 1:20:30
##                    1 h 20 min 30 secs
#@ job_name = mytest
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/mydir
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification=always
#@ notify_user=youremail_at_yoursite.xx
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
#setup of environment
module unload mpi.ibm
module load mpi.intel
export OMP_NUM_THREADS=2
#optional: 
#module load mpi_pinning/hybrid_blocked
mpiexec -n 8000./myprog.exe
