#!/bin/bash
#bash script for strong scaling scenario
#submitcommand="llsubmit"
executable="MarDyn_5683.PAR_RELEASE_AVX-gcc-7-ibmmpi-RMM-DOUBLE"
working_directory="."
processes_per_node=1
time_limit=00:30:00
#should be changed:
#GENERIC_ARG_NAME
#GENERIC_ARG_CLASS
#GENERIC_ARG_TIME
#GENERIC_ARG_NODES
#GENERIC_ARG_TOTAL_TASKS
#GENERIC_ARG_ISLAND_COUNT
#GENERIC_ARG_EXECUTABLE

input_file_name1node="ljfluid1node.xml"
onenodesize=`cat $input_file_name1node | grep "<lx>" | sed "s/>/</g" | cut -d '<' -f 3`

ARG_TIME=$time_limit
ARG_EXECUTABLE=$executable

number_of_nodes="1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 9216"
for i in $number_of_nodes
do

	ARG_NAME="mar-weak-$i-sp-35"
	if (("$i" <= 32))
	then
		ARG_CLASS="test"
		ARG_ISLAND_COUNT=1,1
	fi
	if (($i > 32)) && (($i <= 512))
	then
		ARG_CLASS="general"
		ARG_ISLAND_COUNT=1,1
	fi
	if (($i > 512)) && (($i <= 4096))
	then
		ARG_CLASS="large"
		count=$((i / 512))
		ARG_ISLAND_COUNT=$count,$count
	fi
	if (($i > 4096)) && (($i <= 9216))
	then
		ARG_CLASS="special"
		count=$((i / 512))
		ARG_ISLAND_COUNT=$count,$count
	fi
	ARG_NODES=$i
	ARG_TOTAL_TASKS=$i
	number_of_nodes=$i
	size=`echo "" | awk "END {print $onenodesize * $i ^ (1/3) }"`

	#change size of new input:
	configname="$i-nodes.xml"
	sed -e "s/$onenodesize/$size/g" $input_file_name1node > $configname
	
	#sed everything
	#s/GENERIC_ARG_NAME/$ARG_NAME/g;
	#s/GENERIC_ARG_NAME/$ARG_NAME/g; 
	#s/GENERIC_ARG_CLASS/$ARG_CLASS/g; 
	#s/GENERIC_ARG_TIME/$ARG_TIME/g; 
	#s/GENERIC_ARG_NODES/$ARG_NODES/g; 
	#s/GENERIC_ARG_TOTAL_TASKS/$ARG_TOTAL_TASKS/g; 
	#s/GENERIC_ARG_ISLAND_COUNT/$ARG_ISLAND_COUNT/g; 
	#s/GENERIC_ARG_EXECUTABLE/$ARG_EXECUTABLE/g; 
	
	sed -e "s/GENERIC_ARG_NAME/$ARG_NAME/g; s/GENERIC_ARG_NAME/$ARG_NAME/g; s/GENERIC_ARG_CLASS/$ARG_CLASS/g; s/GENERIC_ARG_TIME/$ARG_TIME/g; s/GENERIC_ARG_NODES/$ARG_NODES/g; s/GENERIC_ARG_TOTAL_TASKS/$ARG_TOTAL_TASKS/g; s/GENERIC_ARG_ISLAND_COUNT/$ARG_ISLAND_COUNT/g; s/GENERIC_ARG_EXECUTABLE/$ARG_EXECUTABLE/g;" generic-script-smuc1.sh > job-$number_of_nodes.sh
#	$submitcommand job$number_of_nodes.sh
done


