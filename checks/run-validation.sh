#!/bin/bash
# Testscript running validation tests according to GitHub Actions but without comparison to master
#
# Copyright (c) 2023      Simon Homes
# Copyright (c) 2024      Christoph Niethammer <niethammer@hlrs.de>
#

repoPath="$PWD/.."  # Path to root directory

export OMP_NUM_THREADS=2

inputlist="${repoPath}/examples/example-list.txt"
logfileName="${repoPath}/examples/run-validation.log"
rm -f ${logfileName}

cd ${repoPath}/examples
./run-examples.sh -v --inputlist "${inputlist}" --logfile "${logfileName}"


