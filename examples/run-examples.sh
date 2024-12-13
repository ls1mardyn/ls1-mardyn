#!/bin/bash
# Testscript running all examples from example-list.txt
#
# Copyright (c) 2017-2024 Christoph Niethammer <niethammer@hlrs.de>
#

# Some default values
MARDYN_EXE="$PWD/../src/MarDyn"
MARDYN_ARGS="--steps 20 -v"
MPIRUN_EXE="mpirun"
MPIRUN_ARGS="-n 4 --oversubscribe"
EXAMPLE_LIST_FILE=${EXAMPLE_LIST_FILE:=example-list.txt}
LOGFILE=${LOGFILE:=$PWD/run-examples.log}

TEMP=$(getopt -o vh --long "help,inputlist:,logfile:,mardyn_args:,mardyn_exe:,mpirun_args:,mpirun_exe:,verbose" -- $@)

if [ $? -ne 0 ]; then
  echo "Error parsing commandline"
  exit 1
fi

eval set -- "$TEMP"
unset TEMP

function print_help() {
cat <<EOF
$(basename $0) is a script to execute MarDyn examples from a inputfile list file

A inputlist file contains one path per line to a configuration file to be passed to MarDyn.
MarDyn is executed with one configuration file at a time.
In the inputlist file empty lines and lines starting with '#' are ignored.

Usage:
   $0 [OPTIONS ...]
   $0 [-h|--help]

Options:
 -i,--inputlist   <FILE>    path to an inputfile list
    --logfile     <FILE>    path for the logfile
    --mpirun_exe  <mpirun>  mpirun command or path to mpirun executable
    --mpirun_args <ARGS>    arguments to be passed to the mpirun command
    --mardyn_exe  <MarDyn>  path to a MarDyn executable
    --mardyn_args <ARGS>    arguments to be passed to MarDyn for each input file
 -v,--verbose               verbose output
 -h,--help                  show this help

EOF
}

while true; do
  case "$1" in
    '--logfile')
      LOGFILE="$2"
      shift 2
      continue
      ;;
    '--mardyn_exe')
       MARDYN_EXE="$2"
       shift 2
       continue
       ;;
    '--mardyn_args')
       MARDYN_ARGS="$2"
       shift 2
       continue
       ;;
    '--mpirun_exe')
       MPIRUN_EXE="$2"
       shift 2
       continue
       ;;
    '--mpirun_args')
       MPIRUN_ARGS="$2"
       shift 2
       continue
       ;;
    '--inputlist')
       EXAMPLE_LIST_FILE="$2"
       shift 2
       continue
       ;;
    '-v'|'--verbose')
       VERBOSE=true
       shift
       continue
       ;;
    '-h'|'--help')
       PRINT_HELP_AND_EXIT=true
       shift
       break
       ;;
   '--')
       shift
       break
       ;;
    * )
       echo "ERROR Unknown Option $1"
       exit 1
       ;;
  esac
done

if [ $PRINT_HELP_AND_EXIT ]; then
    print_help
    exit 0
fi


all_examples=()
failed_examples=()

while IFS= read -r line ; do
  # Skip comment lines, starting with `#`.
  [[ "$line" =~ ^#.*$ ]] && continue
  # Skip empty lines
  [[ "$line" =~ ^$ ]] && continue

  example="$line"
  if [ $VERBOSE ]; then
    echo "Going to run example file $example"
  fi
  all_examples+=("$example")
done < "$EXAMPLE_LIST_FILE"

if [ $VERBOSE ]; then
   echo "Writing test outputs to $LOGFILE"
fi

# definitions for easy color output if we are running in a terminal
if [ -t 0 -a -t 1 -a -t 2 ]; then
  Color_Off='\e[0m'       # Text Reset
  IGreen='\e[0;92m'       # Intense Green
  IRed='\e[0;91m'         # Intense Red
  IMagenta='\e[0;95m'     # magenta
fi

logfile=$LOGFILE
date > $logfile
for example in ${all_examples[@]} ; do
  example_dir=$(dirname $example)
  example_inp_file=$(basename $example)
  echo -e -n "Running example ${IMagenta}$example_inp_file${Color_Off} in ${IMagenta}$example_dir${Color_Off} ... "
  echo "Running example $example_inp_file in $example_dir" >> $logfile
  pushd "$example_dir" >/dev/null
  cmd="$MPIRUN_EXE $MPIRUN_ARGS $MARDYN_EXE $MARDYN_ARGS $example_inp_file"
  echo $cmd >>$logfile
  { $cmd ; } >>$logfile 2>&1
  ret=$?
  if [ $ret -eq 0 ]; then
    echo -e "${IGreen}success${Color_Off}"
    echo "Running example $example_inp_file in $example_dir ... Result: success" >> $logfile
  else
    failed_examples=(${failed_examples[@]} $example)
    echo -e "${IRed}failed${Color_Off}"
    echo "Running example $example_inp_file in $example_dir ... Result: failed" >> $logfile
  fi
  popd >/dev/null
  sleep 1 # sleep for 1 second to make cancelation with CTR-C easier
done

echo "----------------------------------------"
echo "Summary of tested examples:"
echo "----------------------------------------"
echo "num tested examples: ${#all_examples[@]}"
echo "num failed examples: ${#failed_examples[@]}"
echo "----------------------------------------"

date >> $logfile
exit ${#failed_examples[@]}

