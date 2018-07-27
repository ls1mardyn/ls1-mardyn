#!/bin/bash
# Step 3: generate jobscripts hybrid MPI x OpenMP variation of number of MPI ranks
echo "Script $0 generating input."

for ((iMPI=1; iMPI <= 32; iMPI *= 2)) ; 
do
	ARG_MPI=$iMPI
	ARG_OMP=$((32 / iMPI))

	outfile="job-$ARG_MPI-$ARG_OMP.sh"
	cp generic-script.sh $outfile

	#sed everything
	sed -i -e "s/GENERIC_ARG_MPI/$ARG_MPI/g;" $outfile
	sed -i -e "s/GENERIC_ARG_OMP/$ARG_OMP/g;" $outfile
done
