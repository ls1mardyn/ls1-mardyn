#!/bin/bash
numNodes=$1
# Step 3: generate jobscripts hybrid MPI x OpenMP variation of number of MPI ranks
echo "script $0 generating input file for $1 nodes."

time_limit=12:00:00
ARG_TIME=$time_limit
numNodes=$1
ARG_NODES=$numNodes

numMPIs="1 2 4 6 8 12 24 48"
for iMPI in $numMPIs ; 
do
	if (("$iMPI" == 48)); then
		ARG_EXECUTABLE="MarDyn_5713.PAR_RELEASE_AVX-gcc-7-RMM-SINGLE-no-openmp"
	else
		ARG_EXECUTABLE="MarDyn_5713.PAR_RELEASE_AVX-gcc-7-RMM-SINGLE"
	fi
	ARG_TOTAL_TASKS=$((numNodes * iMPI))
	ARG_NAME="mar-h-$numNodes-$iMPI-sp-35"
	ARG_NUMOMP=$((48 / iMPI))
	ARG_PPNODE=$iMPI

	outfilename=job-$numNodes-$iMPI-$ARG_NUMOMP.sh
	cp generic-script-hhen.sh $outfilename

	#sed everything
	sed -i -e "s/GENERIC_ARG_NAME/$ARG_NAME/g;" $outfilename
	sed -i -e "s/GENERIC_ARG_CLASS/$ARG_CLASS/g;" $outfilename
	sed -i -e "s/GENERIC_ARG_TIME/$ARG_TIME/g;" $outfilename
	sed -i -e "s/GENERIC_ARG_NODES/$ARG_NODES/g;" $outfilename
	sed -i -e "s/GENERIC_ARG_TOTAL_TASKS/$ARG_TOTAL_TASKS/g;" $outfilename
	sed -i -e "s/GENERIC_ARG_ISLAND_COUNT/$ARG_ISLAND_COUNT/g;" $outfilename
	sed -i -e "s/GENERIC_ARG_EXECUTABLE/$ARG_EXECUTABLE/g;" $outfilename
	sed -i -e "s/GENERIC_ARG_NUMOMP/$ARG_NUMOMP/g;" $outfilename
	sed -i -e "s/GENERIC_ARG_PPNODE/$ARG_PPNODE/g;" $outfilename

done
