#!/bin/bash
numNodes=$1
# Step 3: generate jobscripts hybrid MPI x OpenMP variation of number of MPI ranks
echo "script $0 generating input file for $1 nodes."

time_limit=05:00:00
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

numMPIs="1 2 4 6 8 12 24 48"
for iMPI in $numMPIs ; 
do
	if (("$iMPI" == 48)); then
		ARG_EXECUTABLE="MarDyn_5713.PAR_RELEASE_AVX-gcc-7-ibmmpi-RMM-SINGLE-no-openmp"
	else
		ARG_EXECUTABLE="MarDyn_5713.PAR_RELEASE_AVX-gcc-7-ibmmpi-RMM-SINGLE"
	fi
	ARG_TOTAL_TASKS=$((numNodes * iMPI))
	ARG_NAME="mar-h-$numNodes-$iMPI-sp-35"
	ARG_NUMOMP=$((48 / iMPI))

	#sed everything
	#s/GENERIC_ARG_NAME/$ARG_NAME/g; 
	#s/GENERIC_ARG_CLASS/$ARG_CLASS/g; 
	#s/GENERIC_ARG_TIME/$ARG_TIME/g; 
	#s/GENERIC_ARG_NODES/$ARG_NODES/g; 
	#s/GENERIC_ARG_TOTAL_TASKS/$ARG_TOTAL_TASKS/g; 
	#s/GENERIC_ARG_ISLAND_COUNT/$ARG_ISLAND_COUNT/g; 
	#s/GENERIC_ARG_EXECUTABLE/$ARG_EXECUTABLE/g; 
	#s/GENERIC_ARG_NUMOMP/$ARG_NUMOMP/g; 

	sed -e "s/GENERIC_ARG_NAME/$ARG_NAME/g; s/GENERIC_ARG_CLASS/$ARG_CLASS/g; s/GENERIC_ARG_TIME/$ARG_TIME/g; s/GENERIC_ARG_NODES/$ARG_NODES/g; s/GENERIC_ARG_TOTAL_TASKS/$ARG_TOTAL_TASKS/g; s/GENERIC_ARG_ISLAND_COUNT/$ARG_ISLAND_COUNT/g; s/GENERIC_ARG_EXECUTABLE/$ARG_EXECUTABLE/g; s/GENERIC_ARG_NUMOMP/$ARG_NUMOMP/g; " generic-script-smuc1.sh > job-$numNodes-$iMPI-$ARG_NUMOMP.sh
done
