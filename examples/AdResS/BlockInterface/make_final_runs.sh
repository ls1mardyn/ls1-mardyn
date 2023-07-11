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
  for n in 2 3
    do
      for type in "${TYPES[@]}"
      do
        for f in "$BUILD_PATH"/src/CP_BI_1024k"${n}"_"${type}"_*1.restart.*(N)
        do
          cp "$f" ./1024k"${n}"
        done
      done
    done
fi
# endregion

# region run
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
for n in 2 3
do
  for type in "${TYPES[@]}"
  do
    rm -f CI_RUN_1024k${n}_"${type}".cmd
    touch CI_RUN_1024k${n}_"${type}".cmd
    echo "#!/bin/bash
#SBATCH -J CI_RUN_1024k${n}_${type}
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

mpiexec -n 8 ./MarDyn ../../examples/AdResS/BlockInterface/1024k${n}/config_${type}_run.xml --loop-abort-time=2370 --final-checkpoint=0
" > CI_RUN_1024k${n}_"${type}".cmd
  done
done
#endregion

#region bench base
for node in 1 4 8 12 16 20 24 28 32 48 64 80 96 112 128
do
  PART="medium"
  if [ $node -lt 7 ]; then PART="small"; fi
  for n in 1 4 16 64 256 1024
  do
    for type in "${TYPES[@]}"
    do
      rm -f CI_BENCH_${n}k_"${type}"_N"${node}".cmd
      touch CI_BENCH_${n}k_"${type}"_N"${node}".cmd
      echo "#!/bin/bash
#SBATCH -J CI_BENCH_${n}k_${type}_N${node}
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

mpiexec -n $((node*4)) ./MarDyn ../../examples/AdResS/BlockInterface/${n}k/config_${type}_run.xml --loop-abort-time=570 --steps=10000 --final-checkpoint=0
" > CI_BENCH_${n}k_"${type}"_N"${node}".cmd
    done
  done
  for n in 2 3
    do
      for type in "${TYPES[@]}"
      do
        rm -f CI_BENCH_1024k${n}_"${type}"_N"${node}".cmd
        touch CI_BENCH_1024k${n}_"${type}"_N"${node}".cmd
        echo "#!/bin/bash
#SBATCH -J CI_BENCH_1024k${n}_${type}_N${node}
#SBATCH -o ../../../log/%x.out
#SBATCH -D ../../../build/src
#SBATCH --partition=$PART
#SBATCH --get-user-env
#SBATCH --mail-type=all
#SBATCH --mem=800mb
#SBATCH --mail-user=hocksa@hsu-hh.de
#SBATCH --export=NONE
#SBATCH --time=00:20:00
#SBATCH --nodes=${node}
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18
#SBATCH --error=../../../log/%x.err

export OMP_NUM_THREADS=18
ulimit -s 1000000

module load slurm_setup
module load mpi/2021.6.0

mpiexec -n $((node*4)) ./MarDyn ../../examples/AdResS/BlockInterface/1024k${n}/config_${type}_run.xml --loop-abort-time=1170 --steps=10000 --final-checkpoint=0
" > CI_BENCH_1024k${n}_"${type}"_N"${node}".cmd
      done
    done
done
#endregion

#region change hybrid
for n in 256 1024
do
  for type in "${TYPES[@]}"
  do
    for percent in {10..90..10}
    do
      rm -f ./${n}k/config_"${type}"_H"${percent}".xml
      cp ./${n}k/config_"${type}"_run.xml ./${n}k/config_"${type}"_H"${percent}".xml
      SIZE=$(sed '27q;d' ./${n}k/config_"${type}"_H"${percent}".xml | sed -r 's#.*<.*>([0-9]*).*</lx>#\1#g' -)
      sed -Ei "s#<steps>100000#<steps>10000#g" ./${n}k/config_"${type}"_H"${percent}".xml
      sed -Ei "s#<writefrequency>100000#<writefrequency>10000#g" ./${n}k/config_"${type}"_H"${percent}".xml
      sed -Ei "s#<outputprefix>BlockInterface_${n}k_${type}#<outputprefix>BlockInterface_${n}k_${type}_H${percent}#g" ./${n}k/config_"${type}"_H"${percent}".xml
      MAX_LEN=$(((SIZE/2)-20))
      LEN=$(((percent/100.0)*MAX_LEN))
      sed -Ei "s#<hybridDimX>.*</hybridDimX>#<hybridDimX>${LEN}</hybridDimX>#g" ./${n}k/config_"${type}"_H"${percent}".xml
    done
  done
done
for n in 2 3
do
  for type in "${TYPES[@]}"
  do
    for percent in {10..90..10}
    do
      rm -f ./1024k${n}/config_"${type}"_H"${percent}".xml
      cp ./1024k${n}/config_"${type}"_run.xml ./1024k${n}/config_"${type}"_H"${percent}".xml
      SIZE=$(sed '27q;d' ./1024k${n}/config_"${type}"_H"${percent}".xml | sed -r 's#.*<.*>([0-9]*).*</lx>#\1#g' -)
      sed -Ei "s#<steps>100000#<steps>10000#g" ./1024k${n}/config_"${type}"_H"${percent}".xml
      sed -Ei "s#<writefrequency>100000#<writefrequency>10000#g" ./1024k${n}/config_"${type}"_H"${percent}".xml
      sed -Ei "s#<outputprefix>BlockInterface_1024k${n}_${type}#<outputprefix>BlockInterface_1024k${n}_${type}_H${percent}#g" ./1024k${n}/config_"${type}"_H"${percent}".xml
      MAX_LEN=$(((SIZE/2)-20))
      LEN=$(((percent/100.0)*MAX_LEN))
      sed -Ei "s#<hybridDimX>.*</hybridDimX>#<hybridDimX>${LEN}</hybridDimX>#g" ./1024k${n}/config_"${type}"_H"${percent}".xml
    done
  done
done


for node in 1 4 8 12 16 20 24 28 32 48 64 80 96 112 128
do
  TIME="10"
  PART="medium"
  LIMIT="570"
  if [ $node -lt 7 ]; then PART="small"; fi
  if [ $node -gt 16 ]; then TIME="05"; LIMIT="270"; fi

  for n in 256 1024
  do
    for type in "${TYPES[@]}"
    do
      for percent in {10..90..10}
      do
        rm -f CI_BENCH_${n}k_"${type}"_H"${percent}"_N"${node}".cmd
        touch CI_BENCH_${n}k_"${type}"_H"${percent}"_N"${node}".cmd
        echo "#!/bin/bash
#SBATCH -J CI_BENCH_${n}k_${type}_H${percent}_N${node}
#SBATCH -o ../../../log/%x.out
#SBATCH -D ../../../build/src
#SBATCH --partition=$PART
#SBATCH --get-user-env
#SBATCH --mail-type=all
#SBATCH --mem=800mb
#SBATCH --mail-user=hocksa@hsu-hh.de
#SBATCH --export=NONE
#SBATCH --time=00:$TIME:00
#SBATCH --nodes=${node}
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18
#SBATCH --error=../../../log/%x.err

export OMP_NUM_THREADS=18
ulimit -s 1000000

module load slurm_setup
module load mpi/2021.6.0

mpiexec -n $((node*4)) ./MarDyn ../../examples/AdResS/BlockInterface/${n}k/config_${type}_H${percent}.xml --loop-abort-time=${LIMIT} --steps=10000 --final-checkpoint=0
" > CI_BENCH_${n}k_"${type}"_H"${percent}"_N"${node}".cmd
      done
    done
  done
  for n in 2 3
    do
      for type in "${TYPES[@]}"
      do
        for percent in {10..90..10}
        do
          rm -f CI_BENCH_1024k${n}_"${type}"_H"${percent}"_N"${node}".cmd
          touch CI_BENCH_1024k${n}_"${type}"_H"${percent}"_N"${node}".cmd
          echo "#!/bin/bash
  #SBATCH -J CI_BENCH_1024k${n}_${type}_H${percent}_N${node}
  #SBATCH -o ../../../log/%x.out
  #SBATCH -D ../../../build/src
  #SBATCH --partition=$PART
  #SBATCH --get-user-env
  #SBATCH --mail-type=all
  #SBATCH --mem=800mb
  #SBATCH --mail-user=hocksa@hsu-hh.de
  #SBATCH --export=NONE
  #SBATCH --time=00:$TIME:00
  #SBATCH --nodes=${node}
  #SBATCH --ntasks-per-node=4
  #SBATCH --cpus-per-task=18
  #SBATCH --error=../../../log/%x.err

  export OMP_NUM_THREADS=18
  ulimit -s 1000000

  module load slurm_setup
  module load mpi/2021.6.0

  mpiexec -n $((node*4)) ./MarDyn ../../examples/AdResS/BlockInterface/1024k${n}/config_${type}_H${percent}.xml --loop-abort-time=${LIMIT} --steps=10000 --final-checkpoint=0
  " > CI_BENCH_1024k${n}_"${type}"_H"${percent}"_N"${node}".cmd
        done
      done
    done
done
# endregion