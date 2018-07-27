#!/bin/bash
#bash script for weak scaling scenario
submitcommand="llsubmit"
executable="MarDyn-gcc72-ibmmpi-14-DOUBLE-5676"
working_directory="."
processes_per_node=2
time_limit=00:45:00
#should be changed:
#GENERIC_ARG_NAME
#GENERIC_ARG_CLASS
#GENERIC_ARG_TIME
#GENERIC_ARG_NODES
#GENERIC_ARG_TOTAL_TASKS
#GENERIC_ARG_ISLAND_COUNT
#GENERIC_ARG_EXECUTABLE

input_file_name1node="config-1-node.xml"
onenodesize=`cat $input_file_name1node | grep "<lx>" | sed "s/>/</g" | cut -d '<' -f 3`
onenodesize_numblocks=`cat $input_file_name1node | grep "<xz>" | sed "s/>/</g" | cut -d '<' -f 3`

echo $onenodesize
#echo $onenodesize_numblocks

ARG_TIME=$time_limit
ARG_EXECUTABLE=$executable

for (( i=1; i<=4096;i=i*2 ))
do
	if (("$i" > 9216))
	then
		i=9216
	fi
	echo $i	
	ARG_NAME="mardyn-weak-$i-bigger-benzen"
	reservation="none"
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
		reservation='srv03-ib\.56\.r'
	fi
	ARG_NODES=$i
#hyperthreading:
	ARG_TOTAL_TASKS=$((2*i))
	number_of_nodes=$i

	numblocks=`echo "" | awk "END {print int($onenodesize_numblocks * $i ^ (1/3)) }"`
	
	size=`echo "" | awk "END {printf \"%.15g\n\", $onenodesize * $numblocks / $onenodesize_numblocks }"`

	#change size of new input:
	configname="$i-nodes.xml"
	
	sed -e "s/$onenodesize/$size/g" -e "s/<xz>$onenodesize_numblocks<\/xz>/<xz>$numblocks<\/xz>/" -e "s/<vapor>$onenodesize_numblocks<\/vapor>/<vapor>$numblocks<\/vapor>/" $input_file_name1node > $configname
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
	if [ "$reservation" == "none" ]
	then
		sed -i.bak "/GENERIC_ARG_RESERVATION/d" job-$number_of_nodes.sh
	else
		sed -i.bak "s/GENERIC_ARG_RESERVATION/$reservation/g" job-$number_of_nodes.sh
	fi
	rm job-$number_of_nodes.sh.bak
#	$submitcommand job-$number_of_nodes.sh
done


