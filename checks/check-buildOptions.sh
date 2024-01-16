#!/bin/bash
# Script to test the building process by altering through different build options
# It uses the default settings and sets the specified options one by one

numJobs=20

checkFolder=$PWD
buildFolder=$PWD/../build_test_$(date '+%Y-%m-%d')

declare -A optionsList

optionsList[ENABLE_MPI]=ON
optionsList[VECTOR_INSTRUCTIONS]=SSE
optionsList[VECTOR_INSTRUCTIONS]=AVX2
optionsList[CMAKE_BUILD_TYPE]=Debug
optionsList[OPENMP]=ON
optionsList[ENABLE_AUTOPAS]=ON
optionsList[ENABLE_UNIT_TESTS]=ON
optionsList[QUICKSCHED]=ON
optionsList[FMM_FFT]=ON
optionsList[WITH_PAPI]=ON
optionsList[ENABLE_ALLLBL]=ON
optionsList[ENABLE_VTK]=ON
optionsList[REDUCED_MEMORY_MODE]=ON


# CMAKE
mkdir -p $buildFolder
cd $buildFolder

echo "Running CMAKE in folder: $buildFolder"

for option in "${!optionsList[@]}"
do
	echo "   Running option: ${option}=${optionsList[${option}]}"
	rm $buildFolder/* -rf
	( CC=mpicc CXX=mpicxx cmake -D${option}=${optionsList[${option}]} .. && make -j${numJobs} ; echo "$?" ) &> $checkFolder/build_test_cmake_${option}.log
done


# MAKE
cd $checkFolder/../src

echo "Running make in src"

for option in "${!optionsList[@]}"
do
	echo "   Running option: ${option}=${optionsList[${option}]}"
	make clean &> /dev/null
	( make -j${numJobs} -f Makefile ${option}=${optionsList[${option}]} ; echo "$?" ) &> $checkFolder/build_test_make_${option}.log
done

