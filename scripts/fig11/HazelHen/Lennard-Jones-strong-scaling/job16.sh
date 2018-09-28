#!/bin/bash
#PBS -N mardyn
#PBS -l nodes=21:ppn=6
#PBS -l walltime=05:00:00
#PBS -o 16nodes.out
#PBS -e 16nodes.err
#PBS -m abe
#PBS -M seckler@in.tum.de,tchipev@in.tum.de
# Change to the directory that the job was submitted from
cd $PBS_O_WORKDIR
# Launch the parallel job to the allocated compute nodes
#printenv | grep PBS
HYPERTHREADING=2
npes=$PBS_NP
ppn=$PBS_NUM_PPN
omps=$((24*$HYPERTHREADING/$ppn))
nodes=$(($npes/$ppn))

targetnodes=16
targetnpes=$(($targetnodes*$ppn))

excludelist="16nodes.out.excludelist"
touch $excludelist
rm $excludelist
touch $excludelist

echo ""
echo "running with $ppn processes per node ($omps omp threads per process) on $nodes nodes."
echo "using $HYPERTHREADING threads per core."
echo ""
export OMP_NUM_THREADS=$omps

#export MPICH_MAX_SHORT_MSG_SIZE=4096
#export MPICH_UNEX_BUFFER_SIZE=8M

module list
#aprun -n $npes -N $ppn -d $OMP_NUM_THREADS -cc depth a.out
executable="MarDyn.PAR_RELEASE_AVX2_OPENMP_5713_gcc"
numExclude=0
# for loop over system sizes
numRepetitions=10
for ((iRepe=0; iRepe <= $numRepetitions; iRepe +=1)) ;
do
    echo "repetition $iRepe"
    schemes="slice c08"
    for iScheme in $schemes ;
    do 
        inputFileName=ljfluid8node-$iScheme.xml
        outputFileName="output-of-jobs/16nodes.out-$iScheme-$iRepe.txt"
        #mpiexec -n $NumProcs ../../$executableName $inputFileName --steps 11 --final-checkpoint=0 >$outputFileName
	#aprun -n $npes -N $ppn -d $OMP_NUM_THREADS -cc depth -j $HYPERTHREADING ../$executable ljfluid8node.xml 50 --final-checkpoint=0
	if (( $numExclude > 0 ))
	then
		aprun --exclude-node-list-file $excludelist -n $targetnpes -N $ppn -d $OMP_NUM_THREADS -cc cpu -j $HYPERTHREADING ../../exec/$executable $inputFileName --steps=11 --final-checkpoint=0 > $outputFileName
	else
		aprun -n $targetnpes -N $ppn -d $OMP_NUM_THREADS -cc cpu -j $HYPERTHREADING ../../exec/$executable $inputFileName --steps=11 --final-checkpoint=0 > $outputFileName
	fi
	

	if cat $outputFileName | grep SlowNode | grep -q "takes longer"
	then
		echo "slow nodes:"
		cat $outputFileName | grep SlowNode | grep "takes longer" | awk -F": Node " '{print $2}' | cut -d':' -f1 | sort -u | sed -e 's/nid0000//g; s/nid000//g; s/nid00//g; s/nid0//g; s/nid//g'
		cat $outputFileName | grep SlowNode | grep "takes longer" | awk -F": Node " '{print $2}' | cut -d':' -f1 | sort -u | sed -e 's/nid0000//g; s/nid000//g; s/nid00//g; s/nid0//g; s/nid//g' >> $excludelist
		numExclude=`cat $excludelist | wc -l`
		if (( $nodes - $numExclude < $targetnodes ))
		then
			echo "too many slow nodes!: $numExclude"
			echo "fixing targetnodes + targetnpes"
			targetnodes=$(($nodes - $numExclude))
			while true
			do
				echo "checking with $targetnodes"
				largestFactor=`factor $targetnodes | rev | cut -d' ' -f1 | rev`
				if (($largestFactor < 20 ))
				then
					echo "works"
					break
				fi
				targetnodes=$(($targetnodes - 1))				
			done
			targetnpes=$(($targetnodes * $ppn))
		fi
		echo "executing with $targetnodes nodes"
	fi
	#aprun -n 48 -N 24 ./exec/MarDyn.PAR_DEBUG_AOS_4251 inp/simple-lj.cfg 10 --final-checkpoint=0
    done
done

