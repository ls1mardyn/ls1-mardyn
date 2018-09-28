#!/bin/bash
#@ job_type = parallel
#@ energy_policy_tag = tchipev00
#@ minimize_time_to_solution = yes
###@ ll_res_id=srv03-ib.56.r
###@ energy_policy_tag = lu78toq.tag_min_t2s
###@ group = vip
###@ set_frequency = 2.3
#@ initialdir=/home/hpc/pr48te/di56giq5/svn/WR_2017/SuperMUC-Phase1/simd-performance
#@ job_name = vec-AVX-1-1-sp-25
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
executableName=MarDyn_5713.SEQ_RELEASE_AVX-gcc-7-nompi-noomp-RMM-SINGLE


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

iVec=AVX
iRc=25
iRepe=2
iSize=32
schemes="slice c08"
echo "scheme: $iScheme"
inputFileName=lj-rc$iRc.xml
outputFileName="output-of-jobs-1-nodes/out-$iVec-$iRc-$iRepe.txt"
./$executableName $inputFileName --steps 9 --final-checkpoint=0 >$outputFileName
