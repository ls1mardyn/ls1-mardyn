#!/usr/bin/zsh

TYPES=('C6H12' 'CH4')
: "${BUILD_PATH:=""}"

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

for n in 1 4 16 64 256 1024
do
  for type in "${TYPES[@]}"
  do
    for f in "$BUILD_PATH"/src/CP_BI_"${n}"k_"${type}"_*1.restart.*(N)
    do
      cp "$f" ./"${n}"k
    done
  done
done

for n in 1 4 16 64 256 1024
do
  for type in "${TYPES[@]}"
  do
    rm -f CI_RUN_${n}k_"${type}".cmd
    touch CI_RUN_${n}k_"${type}".cmd
    echo "#!/bin/bash
#SBATCH -J CI_RUN_${n}k_${type}
#SBATCH -o ../../../log/%x.out
#SBATCH -D ../../../build/src
#SBATCH --partition=medium
#SBATCH --get-user-env
#SBATCH --mail-type=all
#SBATCH --mem=800mb
#SBATCH --mail-user=hocksa@hsu-hh.de
#SBATCH --export=NONE
#SBATCH --time=00:40:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --error=../../../log/%x.err

export OMP_NUM_THREADS=72
ulimit -s 1000000

module load slurm_setup
module load mpi/2021.6.0

mpiexec -n 8 ./MarDyn ../../examples/AdResS/BlockInterface/${n}k/config_${type}_run.xml --loop-abort-time=2370 --final-checkpoint=0
" > CI_RUN_${n}k_"${type}".cmd
  done
done