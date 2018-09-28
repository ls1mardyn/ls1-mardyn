#!/bin/bash
#bash script for strong scaling scenario
submitcommand="llsubmit"
executable="MarDyn-gcc72-ibmmpi-14-DOUBLE-5713-nooverlapping"
working_directory="."
time_limit=1000      #in seconds
#should be changed:
#GENERIC_ARG_NAME
#GENERIC_ARG_CLASS
#GENERIC_ARG_TIME
#GENERIC_ARG_NODES
#GENERIC_ARG_TOTAL_TASKS
#GENERIC_ARG_ISLAND_COUNT
#GENERIC_ARG_EXECUTABLE


ARG_EXECUTABLE=$executable
ARG_TIME=$time_limit

#for (( i=8; i<=16384;i=i*2 ))
for (( i=8; i<=4096;i=i*2 ))
do
	if (("$i" > 9216))
	then
		i=9216
	fi
	echo $i	
	ARG_NAME="mardyn-strong-$i-benzen-c08"
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
	$submitcommand job-$number_of_nodes.sh
done


