#!/bin/bash
#bash script for strong scaling scenario
submitcommand="qsub -W group_list=uhs44130"
executable="MarDyn.PAR_RELEASE_AVX2_OPENMP_5713_gcc"
working_directory="."
processes_per_node=1
hyperthreading=2
time_limit=00:40:00
#should be changed:
#number_of_nodes=1
#output_file_name="${i}nodes.out"
#error_file_name="${i}nodes.err"
#input_file_name="${i}nodes.xml"

input_file_name1node="ljfluid1node.xml"
onenodesize=`cat $input_file_name1node | grep "<lx>" | sed "s/>/</g" | cut -d '<' -f 3`

for (( i=1; i<=9000;i=i*2 ))
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
	input_file_name="${i}nodes.xml"
	output_file_name="${i}nodes.out"
	error_file_name="${i}nodes.err"
	size=`echo "" | awk "END {print $onenodesize * $i ^ (1/3) }"`
	#change size of new input:
	sed -e "s/$onenodesize/$size/g" $input_file_name1node > $input_file_name
	
	
	if (( $number_of_nodes > 4096 ))
	then
		submitcommand="$submitcommand -q XXL"
	fi
	sed -e "s/EXARG/$executable/g; s/WDARG/$working_directory/g; s/IFARG/$input_file_name/g; s/EFARG/$error_file_name/g; s/OFARG/$output_file_name/g; s/NNARG/$number_of_nodes/g; s/TARGETNODES/$i/g; s/PPNARG/$processes_per_node/g; s/HYPARG/$hyperthreading/g; s/TLARG/$time_limit/g" genericJob_tolerant.sh > job$i.sh
	#$submitcommand job$number_of_nodes.sh
done


