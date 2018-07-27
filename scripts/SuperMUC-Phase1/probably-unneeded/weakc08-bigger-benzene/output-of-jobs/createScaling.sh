#!/bin/bash
output="data.txt"
FILES="`ls *.out.* | sort -V`"

echo -e "nodes\tprocesses\topenmp\tgflopmain\tgflopforce" > $output

for f in $FILES
do
	echo ""
	echo "file: $f"
	number_of_nodes=`echo $f | cut -d '.' -f 3`
	number_of_processes=`cat $f | grep "Running with" | grep "processes" | xargs | cut -d ' ' -f 7`
	number_of_openmp=`cat $f | grep "Running with" | grep "OpenMP" | xargs | cut -d ' ' -f 7`
	echo "nodes $number_of_nodes, procs $number_of_processes, openmp $number_of_openmp"
	flop_main=`cat $f | grep "FLOP-rate for main loop" | tail -1 | cut -d ':' -f 2 | cut -d '(' -f 1 | sed -e 's:FLOP/sec::'`
	flop_force=`cat $f | grep "FLOP-rate in force" | tail -1 | cut -d ':' -f 2 | cut -d '(' -f 1 | sed -e 's:FLOP/sec::'`
	

#main
	gflop=`echo $flop_main | cut -d ' ' -f 1`
	unit=`echo $flop_main | cut -d ' ' -f 2`
	echo "$gflop; $unit"
	if [[ $unit == 'T' ]]
	then
		#gflop=`echo "$gflop*1000" | bc`
		gflop=`awk "BEGIN { print $gflop * 1.e3 }"`
	fi
	if [[ $unit == 'P' ]]
	then
		#gflop=`echo "$gflop*1000000" | bc`
		gflop=`awk "BEGIN { print $gflop * 1.e6 }"`
	fi
	gflop_main=$gflop
	
#force
	gflop=`echo $flop_force | cut -d ' ' -f 1`
	unit=`echo $flop_force | cut -d ' ' -f 2`
	echo "$gflop; $unit"
	if [[ $unit == 'T' ]]
	then
		#gflop=`echo "$gflop*1000" | bc`
		gflop=`awk "BEGIN { print $gflop * 1.e3 }"`
	fi
	if [[ $unit == 'P' ]]
	then
		#gflop=`echo "$gflop*1000000" | bc`
		gflop=`awk "BEGIN { print $gflop * 1.e6 }"`
	fi
	gflop_force=$gflop

	echo "gflop_main $gflop_main"
	echo "gflop_force $gflop_force"
	echo -e "$number_of_nodes\t$number_of_processes\t$number_of_openmp\t$gflop_main\t$gflop_force" >> $output
done
