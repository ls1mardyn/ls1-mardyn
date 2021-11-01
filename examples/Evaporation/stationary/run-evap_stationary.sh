#!/bin/bash

MARDYN_EXE="mpirun -np 8 $PWD/../../../src/MarDyn"
MARDYN_OPTIONS="--final-checkpoint=0"
LOGFILE=${LOGFILE:=$PWD/run-evap_stationary.log}


echo "Running example: stationary evaporation" | tee $LOGFILE

cd sim01/run01
echo "Running sim01/run01" | tee -a $LOGFILE
$MARDYN_EXE $MARDYN_OPTIONS config.xml >> $LOGFILE

cd ../../sim01/run02
echo "Running sim01/run02" | tee -a $LOGFILE
$MARDYN_EXE $MARDYN_OPTIONS config.xml >> $LOGFILE

# ---

cd ../../sim02/run01
echo "Running sim02/run01" | tee -a $LOGFILE
$MARDYN_EXE $MARDYN_OPTIONS config.xml >> $LOGFILE

cd ../../sim02/run02
echo "Running sim02/run02" | tee -a $LOGFILE
$MARDYN_EXE $MARDYN_OPTIONS config.xml >> $LOGFILE

cd ../../sim02/run03
echo "Running sim02/run03" | tee -a $LOGFILE
$MARDYN_EXE $MARDYN_OPTIONS config.xml >> $LOGFILE

# ---

cd ../../sim03/run01
echo "Running sim03/run01" | tee -a $LOGFILE
$MARDYN_EXE $MARDYN_OPTIONS config.xml >> $LOGFILE

cd ../../sim03/run02
echo "Running sim03/run02" | tee -a $LOGFILE
$MARDYN_EXE $MARDYN_OPTIONS config.xml >> $LOGFILE

cd ../../sim03/run03
echo "Running sim03/run03" | tee -a $LOGFILE
$MARDYN_EXE $MARDYN_OPTIONS config.xml >> $LOGFILE

echo "Done"

# Clean-up command
#find sim0* -type f -not -name "config.xml" -exec rm {} \;

