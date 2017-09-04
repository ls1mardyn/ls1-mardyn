#!/bin/bash
#bash script for strong scaling scenario
submitcommand="qsub"
executable="MarDyn.PAR_RELEASE_AVX2_OPENMP_5044_gcc_NOGEN_noenergy"
working_directory="."
input_file_name="ljfluid1node.xml"
processes_per_node=2
hyperthreading=1

#should be changed:
#time_limit="00:30:00"
#number_of_nodes=1
#output_file_name="${i}nodes.out"

hours1=4
minutes1=00

for (( i=1; i<=4096;i=i*2 ))
do
	number_of_nodes=$i
	output_file_name="${i}nodes.out"
	error_file_name="${i}nodes.err"
	hours=$((hours1/$i))
	minutes=$((($hours1 % $i*60+$minutes1) / $i))
	if [ $hours -eq 0 ]
	then
		minutes=$(( $minutes < 10 ? 10 : $minutes ))
	fi

	if [ $minutes -lt 10 ]
	then
		minutes="0$minutes"
	fi
	time_limit="0$hours:$minutes:00"
	sed -e "s/EXARG/$executable/g; s/WDARG/$working_directory/g; s/IFARG/$input_file_name/g; s/EFARG/$error_file_name/g; s/OFARG/$output_file_name/g; s/NNARG/$number_of_nodes/g; s/PPNARG/$processes_per_node/g; s/HYPARG/$hyperthreading/g; s/TLARG/$time_limit/g" genericJob.sh > job$number_of_nodes.sh
	$submitcommand job$number_of_nodes.sh
done


