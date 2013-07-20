#!/bin/bash

wget http://www5.in.tum.de/mardyn/scenarios/ljfluid_640k_rc3.tar.gz
tar xfz ljfluid_640k_rc3.tar.gz
rm ljfluid_640k_rc3.tar.gz

# power for processor number
POW=6
START_TIME=$(date +%s)
RETVAL=0

# setup file "current" which will contain the runtime results of this run
# first, print first line with #procs
rm current
echo -n "|| <b> \#procs   " >> current
for ((i = 0; i <= $POW; i=$i+1))
do
    let j=2**$i
    echo -n "| $j" >> current
done
echo "" >> current

# ids of jobs submitted, comma seperated for "squeue"
JOBIDS=""
# ids of jobs submitted, seperated by white space for "scancel"
JOBCANCELIDS=""

# submit batch jobs
for ((i = 0; i <= $POW; i=$i+1))
do  
    let j=2**$i
    JOBID=$(sbatch script-$j.sh | cut -d " " -f 4)
    JOBIDS="$JOBIDS$JOBID,"
    JOBCANCELIDS="$JOBCANCELIDS $JOBID"
done 
echo "waiting for jobs: $JOBIDS"

# wait for completion of all jobs
while [ $(squeue -u lu32reb2 -j $JOBIDS | wc -l) -gt 1 ]
do
    echo "Going to sleep!"
    sleep 500
    # abort if tests take longer than 5 hours
    if [ $(($(date +%s) - $START_TIME)) -gt 18000 ]
    then
	echo "TIMEOUT while executing jobs... TERMINATE!"
        
        for JOB in $JOBCANCELIDS
        do
            $(scancel $JOB)
        done
	exit -1	
    fi
done


# compare runtimes with maximum runtime
# create file with current runtimes
for ((i = 0; i <= $POW; i=$i+1))
do
    let j=2**$i
    sed -n '2p' reference > tmpref
    max=$(cat tmpref | cut -d "|" -f $(($i+4)) )
    ./process-output.py $j.out tmp "" 
    current=$(cat tmp)
    rm tmp
    rm tmpref
    
    # case: output not ok, no compute time found
    if [ "$current" == "0" ];
    then
	echo "NO RUNTIME FOUND!"
        RETVAL=1
    fi

    #case: runtime too high
    if [ "$(echo $max '<' $current | bc -l)" -eq 1 ];
    then
        echo "RUNTIME TOO HIGH! Expected $max but got $current"
    fi
    echo -n "| $current" >> current
done
	
svn commit --username eckhardw-ro --password ra1nING -m "new current runtime for scenario homogeneous-1clj" current 

rm ljfluid_640k_rc3.inp
exit $RETVAL
