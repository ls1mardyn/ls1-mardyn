#!/bin/bash
#bash script to get time values from calculations

output="data.txt"

baseMMUPS="0"
basenodes="0"
echo -e "nodes\titerations\tmolecules\tt-init\tt-main\tt-comp\tt-force\tt-decomp\tgflops-last\tMMUPS\tpareff" > $output
for (( i=1; i<=9096;i=i*2 ))
do
	if (($i > 4096))
	then
		i=7168
	fi
	number_of_nodes=$i
	input_file_name="${i}nodes.out"
	if [ -f $input_file_name ]
	then
		time_init=`cat $input_file_name | grep "Initial IO took" | tail -1 | cut -d ':' -f 3 | sed -e 's/ //g; s/sec//g'`	
		time_main=`cat $input_file_name | grep "Computation in main loop took" | tail -1 | cut -d ':' -f 3 | sed -e 's/ //g; s/sec//g'`
		time_comp=`cat $input_file_name | grep "Computation took" | tail -1 | cut -d ':' -f 3 | sed -e 's/ //g; s/sec//g'`
		time_force=`cat $input_file_name | grep "Force calculation took" | tail -1 | cut -d ':' -f 3 | sed -e 's/ //g; s/sec//g'`
		time_decomp=`cat $input_file_name | grep "Decomposition took" | tail -1 | cut -d ':' -f 3 | sed -e 's/ //g; s/sec//g'`
		gflops_last=`cat $input_file_name | grep "FLOP/sec" | grep main | tail -n 2 | head -n 1 | cut -d':' -f 2 | xargs | cut -d'(' -f 1 | sed -e 's:FLOP/sec::g' | sed -e 's/ //g; s/T/e3/g; s/G//g;s/P/e6/g;s/M/e-3/g'`
		mols=`cat $input_file_name | grep "System contains" | cut -d $'\t' -f 3 | cut -d ' ' -f 3`
		its=`cat $input_file_name  | grep "Simulating " | cut -d $'\t' -f 3  | cut -d ' ' -f 2`
		# million molecule updates per second
		MMUPS=`echo "scale=2 ; $mols * $its / $time_main / 10^6" | bc`

		# calc parallel efficiency
		if [ "$basenodes" -eq "0" ]
		then
			basenodes=$number_of_nodes
			baseMMUPS=$MMUPS	
			pareff="1."
		else
			pareff=`awk "BEGIN {print $MMUPS/$baseMMUPS*$basenodes/$number_of_nodes}"`
			echo $pareff
		fi
	else
		time_init=""
		time_main=""
		time_comp=""
		time_force=""
		time_decomp=""
		gflops_fast=""
		mols=""
		its=""
		gflops_last=""
		MMUPS=""
		pareff=""
	fi

	echo -e "$number_of_nodes\t$its\t$mols\t$time_init\t$time_main\t$time_comp\t$time_force\t$time_decomp\t$gflops_last\t$MMUPS\t$pareff" >> $output
done


