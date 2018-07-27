#!/bin/bash
# job.msub
# job script for Hermit1 (http://www.hlrs.de/systems/platforms/)
# Christoph Niethammer <niethammer@hlrs.de>
#
#PBS -l mppwidth=$NPROC
##PBS -l mppnppn=$NPPN
#PBS -l walltime=1:00:00
#PBS -x
#PBS -N $JOBNAME

#. /opt/modules/default/etc/modules.sh
#module swap PrgEnv-cray PrgEnv-gnu

echo -e "JOBNAME:\t$JOBNAME"
echo -e "NPROC:  \t$NPROC"
#echo -e "NPPN:   \t$NPPN"

INPFILE="M$INP_ID-$INPUT_1R.cfg"

#set -x
cd ${PBS_O_WORKDIR}
pwd
cmd="aprun -n $NPROC ./MarDyn ${INPFILE} 1000"
date +%A,%e.%B.%Y,%H:%M:%S.%N
echo "$cmd"
echo "============================================================"
starttime=`date +"%s"`
$cmd
retval=$?
endtime=`date +"%s"`
echo "============================================================"
date +%A,%e.%B.%Y,%H:%M:%S.%N
difftime_sec=$(( ${endtime}-${starttime} ))
difftime_S=$(( ${difftime_sec}%60 ))
difftime_M=$(( ${difftime_sec}/60 ))
difftime_H=$(( ${difftime_M}/60 ))
difftime_M=$(( ${difftime_M}%60 ))
echo "running for ${difftime_sec} sec = ${difftime_H} h, ${difftime_M} min, ${difftime_S} sec"
if [ ${retval} -gt 0 ]; then
	echo "ERROR (${retval})"
fi

echo "job finished"
