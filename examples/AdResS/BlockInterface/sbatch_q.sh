#!/usr/bin/zsh
: "${INIT_BASE:=0}"
: "${BENCH:=0}"
: "${BENCH_H:=0}"
: "${RUN:=0}"
: "${INIT_XS:=0}"
: "${BENCH_XS:=0}"

TYPES=('C6H12' 'CH4')

if [[ $INIT_BASE -gt 0 ]]
then
  for n in 1 4 16 64 256 1024
  do
    for type in "${TYPES[@]}"
    do
      sbatch ./CI_INIT_${n}k_"${type}"_0.cmd
      sbatch ./CI_INIT_${n}k_"${type}"_1.cmd
    done
  done
fi

if [[ $RUN -gt 0 ]]
then
  for n in 1 4 16 64 256 1024
  do
    for type in "${TYPES[@]}"
    do
      sbatch ./CI_RUN_${n}k_"${type}".cmd
    done
  done
fi

if [[ $BENCH -gt 0 ]]
then
  for node in 1 4 8 12 16 20 24 28 32 48 64 80 96 112 128
  do
    for n in 1 4 16 64 256 1024
    do
      for type in "${TYPES[@]}"
      do
        sbatch ./CI_BENCH_${n}k_"${type}"_N"${node}".cmd
      done
    done
  done
fi

if [[ $BENCH_H -gt 0 ]]
then
  for node in 1 4 8 12 16 20 24 28 32 48 64 80 96 112 128
  do
    for n in 1 4 16 64 256 1024
    do
      for type in "${TYPES[@]}"
      do
        for percent in {10..90..10}
        do
          sbatch ./CI_BENCH_${n}k_"${type}"_H"${percent}"_N"${node}".cmd
        done
      done
    done
  done
fi

if [[ $INIT_XS -gt 0 ]]
then
  for size in {0..1000..100}
  do
    for type in "${TYPES[@]}"
    do
      sbatch ./CI_INIT_XS"${size}"_"${type}"_0.cmd
      sbatch ./CI_INIT_XS"${size}"_"${type}"_1.cmd
    done
  done
fi

if [[ $BENCH_XS -gt 0 ]]
then
  for node in 1 4 8 12 16 20 24 28 32 48 64 80 96 112 128
  do
    for size in {0..1000..100}
    do
      for type in "${TYPES[@]}"
      do
        sbatch ./CI_BENCH_XS"${size}"_"${type}"_N"${node}".cmd
      done
    done
  done
fi