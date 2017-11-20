#!/bin/bash
numNodes=1
ARG_SIZE=32
# Step 3: generate jobscripts hybrid MPI x OpenMP variation of number of MPI ranks
echo "script $0 generating input file for $numNodes nodes for size $ARG_SIZE."

time_limit=00:30:00
ARG_TIME=$time_limit
ARG_NODES=$numNodes

ARG_CLASS="test"
ARG_ISLAND_COUNT=1,1

iMPI=$numNodes
for ((iRepe=0; iRepe <= 2; iRepe +=1)) ; 
do
	ARG_REPE=$iRepe
	iOMP=1
	vecMode="AVX SSE SOA"
	for iVec in $vecMode ; 
	do
		ARG_VEC=$iVec
		ARG_EXECUTABLE="MarDyn_5713.SEQ_RELEASE_$iVec-gcc-7-nompi-noomp-RMM-SINGLE"
		rcMode="25 35 50"
		for iRc in $rcMode ;
		do
			ARG_RC=$iRc
			ARG_TOTAL_TASKS=$((numNodes * iMPI))
			ARG_NAME="vec-$iVec-$numNodes-$iOMP-sp-$iRc"
			ARG_NUMOMP=$iOMP

			outfile="job-$iVec-$iRc-$iRepe.sh"
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
			sed -i -e "s/GENERIC_ARG_RC/$ARG_RC/g;" $outfile
			sed -i -e "s/GENERIC_ARG_VEC/$ARG_VEC/g;" $outfile
		done
	done
done
