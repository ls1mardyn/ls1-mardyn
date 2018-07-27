#!/bin/bash
numNodes=$1
ARG_SIZE=$2
# Step 3: generate jobscripts hybrid MPI x OpenMP variation of number of MPI ranks
echo "script $0 generating input file for $1 nodes for size $2."

time_limit=00:30:00
ARG_TIME=$time_limit
numNodes=$1
ARG_NODES=$numNodes

if (("$numNodes" <= 32))
then
	ARG_CLASS="test"
	ARG_ISLAND_COUNT=1,1
fi
if (($numNodes > 32)) && (($numNodes <= 512))
then
	ARG_CLASS="general"
	ARG_ISLAND_COUNT=1,1
fi
if (($numNodes > 512)) && (($numNodes <= 4096))
then
	ARG_CLASS="large"
	count=$((numNodes / 512))
	ARG_ISLAND_COUNT=$count,$count
fi
if (($numNodes > 4096)) && (($numNodes <= 9216))
then
	ARG_CLASS="special"
	count=$((numNodes / 512))
	ARG_ISLAND_COUNT=$count,$count
fi

iMPI=$numNodes
for ((iRepe=0; iRepe <= 4; iRepe +=1)) ; 
do
	ARG_REPE=$iRepe
	for ((iOMP=1; iOMP <= 32; iOMP *= 2)) ; 
	do
		if (("$iMPI" == 1)); then
			ARG_EXECUTABLE="MarDyn_5713.SEQ_RELEASE_AVX-gcc-7-nompi-RMM-SINGLE"
		else
			ARG_EXECUTABLE="MarDyn_5713.PAR_RELEASE_AVX-gcc-7-ibmmpi-RMM-SINGLE"
		fi
		ARG_TOTAL_TASKS=$((numNodes * iMPI))
		ARG_NAME="omp-strong-$numNodes-$iOMP-sp-35"
		ARG_NUMOMP=$iOMP

		outfile="job-$numNodes-$iMPI-$ARG_NUMOMP-$iRepe.sh"
		cp generic-script-smuc1.sh $outfile

		#sed everything
		sed -i -e "s/GENERIC_ARG_NAME/$ARG_NAME/g;" $outfile
		sed -i -e "s/GENERIC_ARG_CLASS/$ARG_CLASS/g;" $outfile
		sed -i -e "s/GENERIC_ARG_TIME/$ARG_TIME/g;" $outfile
		sed -i -e "s/GENERIC_ARG_NODES/$ARG_NODES/g;" $outfile
		sed -i -e "s/GENERIC_ARG_TOTAL_TASKS/$ARG_TOTAL_TASKS/g;" $outfile
		sed -i -e "s/GENERIC_ARG_ISLAND_COUNT/$ARG_ISLAND_COUNT/g;" $outfile
		sed -i -e "s/GENERIC_ARG_EXECUTABLE/$ARG_EXECUTABLE/g;" $outfile
		sed -i -e "s/GENERIC_ARG_NUMOMP/$ARG_NUMOMP/g;" $outfile
		sed -i -e "s/GENERIC_ARG_SIZE/$ARG_SIZE/g;" $outfile
		sed -i -e "s/GENERIC_ARG_REPE/$ARG_REPE/g;" $outfile
	done
done
