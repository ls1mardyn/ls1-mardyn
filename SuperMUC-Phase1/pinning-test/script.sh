#!/bin/bash
#@ job_type = parallel
#@ energy_policy_tag = tchipev00
#@ minimize_time_to_solution = yes
#@ initialdir=/home/hpc/pr48te/di56giq5/svn/WR_2017/SuperMUC-Phase1/pinning-test
#@ job_name = mar-pinning-test
#@ class = test
#@ node_usage = not_shared
#@ wall_clock_limit = 00:30:00
#@ network.MPI = sn_all,not_shared,us
#@ notification = always
#@ notify_user = tchipev@in.tum.de
#@ output = jobout/LLRUN.out.$(jobid)
#@ error =  jobout/LLRUN.err.$(jobid)
#@ node = 8
#@ island_count=1,1
#@ total_tasks = 8
#@ queue



### @ Other people using this script - change the Email address line! notify_user = xx@in.tum.de


exe="MarDyn_5682.PAR_RELEASE_AVX-gcc-7-ibmmpi"



#-------------------------------------
#Settings for POE
#-------------------------------------
. /etc/profile
. /etc/profile.d/modules.sh
module unload mpi.ibm
module load mpi.ibm/1.4_gcc
export MP_SINGLE_THREAD=no
export OMP_NUM_THREADS=32
# Pinning, this will use the physical and virtual cores
export MP_TASK_AFFINITY=cpu:$OMP_NUM_THREADS

#module load likwid
#likwid-topology >out-likwid-top.txt

#mpiexec -n 8 ./$exe 8max.xml --steps 6 --final-checkpoint=0 > out-default-pinning.txt

#export OMP_PROC_BIND=close
#export OMP_PLACES=threads
#mpiexec -n 8 ./$exe 8max.xml --steps 6 --final-checkpoint=0 > out-close-threads.txt

#export OMP_PROC_BIND=close
#export OMP_PLACES=cores
#mpiexec -n 8 ./$exe 8max.xml --steps 6 --final-checkpoint=0 > out-close-cores.txt

export OMP_PLACES="{0},{16},{1},{17},{2},{18},{3},{19},{4},{20},{5},{21},{6},{22},{7},{23},{8},{24},{9},{25},{10},{26},{11},{27},{12},{28},{13},{29},{14},{30},{15},{31}"
export OMP_PROC_BIND=true
mpiexec -n 8 ./$exe 8max.xml --steps 6 --final-checkpoint=0 > out-manual-true.txt
