#!/bin/sh
#
# build several configurations and run unit tests

# TODO: add the (four) vectorization modes (what about CFGs?)
BUILD1="PARTYPE=SEQ TARGET=DEBUG   UNIT_TESTS=1 VTK=1"
BUILD2="PARTYPE=SEQ TARGET=RELEASE UNIT_TESTS=1"
BUILD3="PARTYPE=SEQ TARGET=RELEASE UNIT_TESTS=0"
BUILD4="PARTYPE=PAR TARGET=DEBUG   UNIT_TESTS=0"
BUILD5="PARTYPE=PAR TARGET=DEBUG   UNIT_TESTS=1"
BUILD6="PARTYPE=PAR TARGET=RELEASE UNIT_TESTS=1 VTK=1"

# list the expected return values for each configuration above (note: make returns 2 if mardyn returns with errors)
RETVAL1=0
RETVAL2=0
RETVAL3=2
RETVAL4=2
RETVAL5=0
RETVAL6=0

NUM_JOBS=3

#LOG_FILE=/dev/null
LOG_FILE=test_build.out
TEST_LOG_FILE=test_build_results.txt

RETURN_VALUE=0

cd src
rm $LOG_FILE
rm $TEST_LOG_FILE
echo "Protocoll of testing build and tests. " > $LOG_FILE

for i in $(seq 6) 
do
  TMP=BUILD$i
  eval ARG=\$$TMP
  
  TMP=RETVAL$i
  eval EXPECTED_RESULT=\$$TMP
 
  echo "===== BUILD $i =====" >> $LOG_FILE  
  (make -s -f ../makefile/Makefile $ARG cleanall) >> $LOG_FILE 2>&1
  (make -s -f ../makefile/Makefile $ARG -j$NUM_JOBS) >> $LOG_FILE 2>&1
  if [ $? -ne 0 ]
    then echo "build   $ARG     FAILED!";
    RETURN_VALUE=1;
  else
    echo "build   $ARG     OK!"
    (make -s -f ../makefile/Makefile $ARG test) >> $TEST_LOG_FILE 2>&1
    eval RETVAL=$?
    if [ $RETVAL -eq $EXPECTED_RESULT ]
      then echo "Tests OK!";
    else
      echo "Tests FAILED!";
      echo Expected Value: $EXPECTED_RESULT;
      echo Return Value: $RETVAL;
      RETURN_VALUE=1;
    fi
  fi 
done

exit $RETURN_VALUE
