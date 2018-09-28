#!/bin/bash
#PBS -W group_list=uhs44130
#PBS -N mar-os-1-6
#PBS -m abe
#PBS -M tchipev@in.tum.de
#PBS -l nodes=1:ppn=6
#PBS -l walltime=1200
#PBS -o output-of-jobs-1-nodes/
#PBS -e output-of-jobs-1-nodes/
#PBS -p 0

# Change to the directory that the job was submitted from
cd $PBS_O_WORKDIR

# Launch the parallel job to the allocated compute nodes
#printenv | grep PBS
omps=6
nodes=1
echo ""
echo "running with $ppn processes per node ($omps omp threads per process) on $nodes nodes."
echo ""
export OMP_NUM_THREADS=$omps

module list

executable="MarDyn_5713.SEQ_RELEASE_AVX-gcc-7-RMM-SINGLE-no-mpi"

# for loop over system sizes
for ((iRepe=0; iRepe <= 4; iRepe +=1)) ; 
do
	echo "repetition $iRepe"
	#for loop over schemes
	schemes="slice c08"
	for iScheme in $schemes ;
	do
		echo "scheme: $iScheme"
		inputFileName=lj-1-64-$iScheme.xml
		outputFileName="output-of-jobs-1-nodes/out-$omps-$iScheme-$iRepe.txt"
		aprun --p-state 2500000 -n 1 -d $OMP_NUM_THREADS -cc cpu -j 1 ../$executable $inputFileName --steps 11 --final-checkpoint=0 > $outputFileName
	done
done
