#!/bin/bash
#bash script for strong scaling scenario
submitcommand="qsub -W group_list=uhs44130"
executable="MarDyn.PAR_RELEASE_AVX2_OPENMP_5713_gcc"
working_directory="."
processes_per_node=6
hyperthreading=2

input_file_name="ljfluid8node.xml"
#should be changed:
#time_limit="00:30:00"
#number_of_nodes=1
#output_file_name="${i}nodes.out"
#error_file_name="${i}nodes.err"
#input_file_name="${i}nodes.xml"

hours1=80
minutes1=00

for (( i=8; i<=9000;i=i*2 ))
do
	puffer=5

	if (( $i == 4096 ))
	then
		i=4096
		puffer=0
	fi

	if (( $i > 4096 ))
	then
		i=7168
		puffer=1
	fi


	number_of_nodes=$(($i+$puffer))
	output_file_name="${i}nodes.out"
	error_file_name="${i}nodes.err"
	hours=$((hours1/$i))
	minutes=$((($hours1 % $i*60+$minutes1) / $i))
	if [ $hours -eq 0 ]
	then
		minutes=$(( $minutes < 30 ? 30 : $minutes ))
	fi

	if [ $minutes -lt 10 ]
	then
		minutes="0$minutes"
	fi
	time_limit="0$hours:$minutes:00"
	if (( $number_of_nodes > 4096 ))
	then
		submitcommand="$submitcommand -q XXL"
	fi
	sed -e "s/EXARG/$executable/g; s/WDARG/$working_directory/g; s/IFARG/$input_file_name/g; s/EFARG/$error_file_name/g; s/OFARG/$output_file_name/g; s/NNARG/$number_of_nodes/g; s/TARGETNODES/$i/g; s/PPNARG/$processes_per_node/g; s/HYPARG/$hyperthreading/g; s/TLARG/$time_limit/g" genericJob_tolerant.sh > job$i.sh
	#$submitcommand job$number_of_nodes.sh
done
