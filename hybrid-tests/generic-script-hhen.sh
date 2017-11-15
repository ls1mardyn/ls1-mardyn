#!/bin/bash
#PBS -N GENERIC_ARG_NAME
#PBS -l nodes=GENERIC_ARG_NODES:ppn=GENERIC_ARG_PPNODE
#PBS -l walltime=GENERIC_ARG_TIME
#PBS -o output-of-jobs-GENERIC_ARG_NODES-nodes/PBS.out.GENERIC_ARG_TOTAL_TASKS.$(PBS_JOBID)
#PBS -e output-of-jobs-GENERIC_ARG_NODES-nodes/PBS.err.GENERIC_ARG_TOTAL_TASKS.$(PBS_JOBID)

# Change to the directory that the job was submitted from
cd $PBS_O_WORKDIR

# Launch the parallel job to the allocated compute nodes
#printenv | grep PBS
npes=$PBS_NP
ppn=$PBS_NUM_PPN
omps=$((48/$ppn))
nodes=$(($npes/$ppn))
echo ""
echo "running with $ppn processes per node ($omps omp threads per process) on $nodes nodes."
echo ""
export OMP_NUM_THREADS=$omps

module list
#aprun -n $npes -N $ppn -d $OMP_NUM_THREADS -cc depth a.out

executable="GENERIC_ARG_EXECUTABLE"

# for loop over system sizes
for ((iRepe=0; iRepe <= 4; iRepe +=1)) ; 
do
	echo "repetition $iRepe"
	for ((iSize=4096; iSize >= 1; iSize /= 2)) ;
	do
		echo "system size: $iSize"
		#for loop over schemes
		schemes="slice c08"
		for iScheme in $schemes ;
		do
			echo "scheme: $iScheme"
			inputFileName=GENERIC_ARG_NODES-nodes-$iSize-$iScheme.xml
			outputFileName="output-of-jobs-GENERIC_ARG_NODES-nodes/out-$NumProcs-$iSize-$iScheme-$iRepe.txt"
			aprun -n $npes -N $ppn -d $OMP_NUM_THREADS -cc depth -j 2 ../$executable $inputFileName 11 --final-checkpoint=0 > $outputFileName
		done
	done
done
