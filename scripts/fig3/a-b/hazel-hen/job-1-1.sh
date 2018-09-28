#!/bin/bash
#PBS -W group_list=uhs44130
#PBS -N mar-simd
#PBS -m abe
#PBS -M tchipev@in.tum.de
#PBS -l nodes=1:ppn=1
#PBS -q test
#PBS -l walltime=1500
#PBS -o output-of-jobs/
#PBS -e output-of-jobs/
#PBS -p 0

# Change to the directory that the job was submitted from
cd $PBS_O_WORKDIR

module list

args="lj-1-64-slice.xml --steps=11 --final-checkpoint=0"

modi="AVX2 AVX SSE SOA"

export OMP_NUM_THREADS=1 
for iVec in $modi ; 
do
	aprun -n 1 -d 1 -cc cpu -j 1 ../MarDyn_5713.SEQ_RELEASE_$iVec $args > output-of-jobs/out-$iVec.txt
done

