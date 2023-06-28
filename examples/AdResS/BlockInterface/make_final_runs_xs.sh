#!/usr/bin/zsh

TYPES=('C6H12' 'CH4')
: "${BUILD_PATH:=""}"
: "${COPY:=0}"

# region path
if [[ BUILD_PATH -eq "" ]]
then
  if [ -d "../../../cmake-build-debug" ]
  then
    BUILD_PATH="../../../cmake-build-debug"
  else
    if [ -d "../../../build" ]
    then
      BUILD_PATH="../../../build"
    else
      exit 255
    fi
  fi
fi
# endregion

# region copy
if [ $COPY -gt 0 ]
then
  for size in {100..1000..100}
  do
    for type in "${TYPES[@]}"
    do
      for f in "$BUILD_PATH"/src/CP_BI_XS"${size}"_"${type}"_*1.restart.*(N)
      do
        cp "$f" ./xScale
      done
    done
  done
fi
# endregion

#region bench base
for node in 1 4 8 12 16 20 24 28 32 48 64 80 96 112 128
do
  PART="medium"
  if [ $node -lt 7 ]; then PART="small"; fi
  for size in {100..1000..100}
  do
    for type in "${TYPES[@]}"
    do
      rm -f CI_BENCH_XS"${size}"_"${type}"_N"${node}".cmd
      touch CI_BENCH_XS"${size}"_"${type}"_N"${node}".cmd
      echo "#!/bin/bash
#SBATCH -J CI_BENCH_XS${size}_${type}_N${node}
#SBATCH -o ../../../log/%x.out
#SBATCH -D ../../../build/src
#SBATCH --partition=$PART
#SBATCH --get-user-env
#SBATCH --mail-type=all
#SBATCH --mem=800mb
#SBATCH --mail-user=hocksa@hsu-hh.de
#SBATCH --export=NONE
#SBATCH --time=00:10:00
#SBATCH --nodes=${node}
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18
#SBATCH --error=../../../log/%x.err

export OMP_NUM_THREADS=18
ulimit -s 1000000

module load slurm_setup
module load mpi/2021.6.0

mpiexec -n $((node*4)) ./MarDyn ../../examples/AdResS/BlockInterface/xScale/config_${type}_XS${size}_run.xml --loop-abort-time=570 --steps=10000 --final-checkpoint=0
" > CI_BENCH_XS"${size}"_"${type}"_N"${node}".cmd
    done
  done
done
#endregion