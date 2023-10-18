#!/bin/bash
# Testscript running validation tests according to GitHub Actions but without comparison to master
#
# Copyright (c) 2023 Simon Homes
#

repoPath="$PWD/.."  # Path to root directory

export OMP_NUM_THREADS=2

logfileName="${repoPath}/examples/run-validation.log"
rm -f ${logfileName}

IFS=$'\n'
for i in $(cat "${repoPath}/examples/example-list.txt" )
do
	# skip if comment or empty line
	if [[ $i == \#* || -z "$i" ]]
	then
		continue
	fi
	cd ${repoPath}/examples/$(dirname $i)

	# run the examples
	printf "Running example: ${i} ... " | tee -a ${logfileName}
	EXE=${repoPath}/build/src/MarDyn

	mpirun --oversubscribe -np 4 ${EXE} $(basename $i) --steps=20 | tee -a ${logfileName}
	printf "done\n"
done

