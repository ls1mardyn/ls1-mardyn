#!/bin/bash
#bash script to get time values from calculations

output="data.txt"

echo -e "#nodes\tt-init\tt-main\tt-comp\tt-force\tt-decomp" > $output
for (( i=1; i<=4096;i=i*2 ))
do
	number_of_nodes=$i
	input_file_name="${i}nodes.out"
	if [ -f $input_file_name ]
	then
		time_init=`cat $input_file_name | grep "Initial IO took" | cut -d ':' -f 3 | sed -e 's/ //g; s/sec//g'`	
		time_main=`cat $input_file_name | grep "Computation in main loop took" | tail -1 | cut -d ':' -f 3 | sed -e 's/ //g; s/sec//g'`
		time_comp=`cat $input_file_name | grep "Computation took" | tail -1 | cut -d ':' -f 3 | sed -e 's/ //g; s/sec//g'`
		time_force=`cat $input_file_name | grep "Force calculation took" | tail -1 | cut -d ':' -f 3 | sed -e 's/ //g; s/sec//g'`
		time_decomp=`cat $input_file_name | grep "Decomposition took" | tail -1 | cut -d ':' -f 3 | sed -e 's/ //g; s/sec//g'`
	else
		time_init=""
		time_main=""
		time_comp=""
		time_force=""
		time_decomp=""
	fi

	echo -e "$number_of_nodes\t$time_init\t$time_main\t$time_comp\t$time_force\t$time_decomp" >> $output
done


