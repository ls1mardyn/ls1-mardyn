#!/bin/bash
#PBS -W group_list=uhs44130
#PBS -N mar-h-1-24
#PBS -m abe
#PBS -M tchipev@in.tum.de
#PBS -l nodes=1:ppn=24
#PBS -l walltime=00:25:00
#PBS -o output-of-jobs-64-nodes/
#PBS -e output-of-jobs-64-nodes/
#PBS -p 60

# Change to the directory that the job was submitted from
cd $PBS_O_WORKDIR

# Launch the parallel job to the allocated compute nodes
#printenv | grep PBS
npes=24
ppn=24
nodes=$(($npes/$ppn))
echo ""
echo "running with $ppn processes per node ($omps omp threads per process) on $nodes nodes."
echo ""
#export OMP_NUM_THREADS=$omps

NumProcs=$npes

module list
#aprun --p-state 2500000 -n $npes -N $ppn -d $OMP_NUM_THREADS -cc depth a.out

executable="MarDyn_5713.PAR_RELEASE_AVX-gcc-7-RMM-SINGLE-no-openmp"

# for loop over system sizes
for ((iRepe=0; iRepe <= 4; iRepe +=1)) ; 
do
	echo "repetition $iRepe"
	#for ((iSize=524288; iSize >= 1; iSize /= 2)) ;
	#do
		iSize=64
		echo "system size: $iSize"
		#for loop over schemes
		schemes="slice"
		for iScheme in $schemes ;
		do
			echo "scheme: $iScheme"
			inputFileName=lj-1-$iSize-$iScheme.xml
			outputFileName="output-of-jobs-64-nodes/out-$NumProcs-$iSize-$iScheme-$iRepe.txt"
			aprun --p-state 2500000 -n $npes -N $ppn -cc cpu -j 1 ../$executable $inputFileName --steps 11 --final-checkpoint=0 > $outputFileName
		done
	#done
done
