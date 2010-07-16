#!/bin/sh
#
# build several configurations and run unit tests

BUILD1="PARTYPE=SEQ TARGET=DEBUG TESTS=1"
BUILD2="PARTYPE=SEQ TARGET=RELEASE TESTS=1"
BUILD3="PARTYPE=SEQ TARGET=RELEASE TESTS=0"
BUILD4="PARTYPE=PAR TARGET=DEBUG TESTS=0"
BUILD5="PARTYPE=PAR TARGET=DEBUG TESTS=1"
BUILD6="PARTYPE=PAR TARGET=RELEASE TESTS=1"

NUM_JOBS=3

cd src

for i in $(seq 6) 
do
  TMP=BUILD$i
  eval ARG=\$$TMP
  (make -s -f ../makefile/Makefile $ARG clean) > /dev/null 2>&1
  (make -s -f ../makefile/Makefile $ARG -j$NUM_JOBS) > /dev/null 2>&1
  if [ $? -ne 0 ]
    then echo "build   $ARG     FAILED!";
  else
    echo "build   $ARG     OK!"
    (make -s -f ../makefile/Makefile $ARG test) > /dev/null 2>&1
    if [ $? -eq 0 ]
      then echo "Tests OK!";
    else
      echo "Tests FAILED!";
    fi
  fi 
done
