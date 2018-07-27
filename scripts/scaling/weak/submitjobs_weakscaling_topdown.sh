#!/bin/bash
#bash script for strong scaling scenario
submitcommand="qsub"
executable="MarDyn.PAR_RELEASE_AVX2_OPENMP_5044_gcc_NOGEN_noenergy"
working_directory="."
processes_per_node=2
hyperthreading=1
time_limit=01:00:00
#should be changed:
#number_of_nodes=1
#output_file_name="${i}nodes.out"
#error_file_name="${i}nodes.err"
#input_file_name="${i}nodes.xml"

input_file_namemaxnode="ljfluid7168node.xml"
maxnodes=7168
maxnodesize=`cat $input_file_namemaxnode | grep "<lx>" | sed "s/>/</g" | cut -d '<' -f 3`
echo "maxnodesize ($maxnodes nodes)=$maxnodesize"
for (( i=1; i<=4096;i=i*2 ))
do
	number_of_nodes=$i
	input_file_name="${i}nodes.xml"
	output_file_name="${i}nodes.out"
	error_file_name="${i}nodes.err"
	size=`echo "" | awk "END {print $maxnodesize * ($i/$maxnodes) ^ (1/3)} "`
	#change size of new input:
	sed -e "s/$maxnodesize/$size/g" $input_file_namemaxnode > $input_file_name
	
	sed -e "s/EXARG/$executable/g; s/WDARG/$working_directory/g; s/IFARG/$input_file_name/g; s/EFARG/$error_file_name/g; s/OFARG/$output_file_name/g; s/NNARG/$number_of_nodes/g; s/PPNARG/$processes_per_node/g; s/HYPARG/$hyperthreading/g; s/TLARG/$time_limit/g" genericJob.sh > job$number_of_nodes.sh
	$submitcommand job$number_of_nodes.sh
done

number_of_nodes=$maxnodes
input_file_name="${maxnodes}nodes.xml"
output_file_name="${maxnodes}nodes.out"
error_file_name="${maxnodes}nodes.err"
size=$maxnodesize
#change size of new input:
sed -e "s/$maxnodesize/$size/g" $input_file_namemaxnode > $input_file_name

sed -e "s/EXARG/$executable/g; s/WDARG/$working_directory/g; s/IFARG/$input_file_name/g; s/EFARG/$error_file_name/g; s/OFARG/$output_file_name/g; s/NNARG/$number_of_nodes/g; s/PPNARG/$processes_per_node/g; s/HYPARG/$hyperthreading/g; s/TLARG/$time_limit/g" genericJob.sh > job$number_of_nodes.sh
$submitcommand -q XXL job$number_of_nodes.sh
