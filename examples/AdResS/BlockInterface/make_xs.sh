#!/usr/bin/zsh

: "${CLEAN:=1}"
TYPES=('C6H12' 'CH4')

if [[ $CLEAN -gt 0 ]]
then
  for size in {0..1000..100}
  do
    for type in "${TYPES[@]}"
    do
      rm -f ./xScale/config_"${type}"_XS"${size}"_0.xml
      rm -f ./xScale/config_"${type}"_XS"${size}"_1.xml
      rm -f ./CI_INIT_XS"${size}"_"${type}"_0.cmd
      rm -f ./CI_INIT_XS"${size}"_"${type}"_1.cmd
      rm -f ./xScale/config_"${type}"_XS"${size}"_run.xml
    done
  done
  exit 0
fi


for size in {0..1000..100}
do
  for type in "${TYPES[@]}"
  do
    rm -f ./xScale/config_"${type}"_XS"${size}"_0.xml
    rm -f ./xScale/config_"${type}"_XS"${size}"_1.xml
    cp ./xScale/config_"${type}".xml ./xScale/config_"${type}"_XS"${size}"_0.xml
    cp ./xScale/config_"${type}".xml ./xScale/config_"${type}"_XS"${size}"_1.xml

    MSIZE=$(sed '27q;d' ./xScale/config_"${type}".xml | sed -r 's#.*<.*>([0-9]*).*</lx>#\1#g' -)
    sed -i -r "s#<lx unit=\"reduced\">.*</lx>#<lx unit=\"reduced\">$((MSIZE/2+size))\.</lx>#g" ./xScale/config_"${type}"_XS"${size}"_0.xml
    sed -i -r "s#<lx unit=\"reduced\">.*</lx>#<lx unit=\"reduced\">$((MSIZE/2+size))\.</lx>#g" ./xScale/config_"${type}"_XS"${size}"_1.xml
    sed -i -n '60,81!p' ./xScale/config_"${type}"_XS"${size}"_0.xml
    sed -i -n '38,59!p' ./xScale/config_"${type}"_XS"${size}"_1.xml
    sed -i -n '80,84!p' ./xScale/config_"${type}"_XS"${size}"_0.xml
    sed -i -n '80,84!p' ./xScale/config_"${type}"_XS"${size}"_1.xml
    sed -i -n '87,96!p' ./xScale/config_"${type}"_XS"${size}"_0.xml
    sed -i -n '87,96!p' ./xScale/config_"${type}"_XS"${size}"_1.xml
    sed -i -r "s#<lower> <x>.*</x>#<lower> <x>0</x>#g" ./xScale/config_"${type}"_XS"${size}"_0.xml
    sed -i -r "s#<upper> <x>.*</x>#<upper> <x>$((MSIZE/2+size))</x>#g" ./xScale/config_"${type}"_XS"${size}"_0.xml
    sed -i -r "s#<lower> <x>.*</x>#<lower> <x>0</x>#g" ./xScale/config_"${type}"_XS"${size}"_1.xml
    sed -i -r "s#<upper> <x>.*</x>#<upper> <x>$((MSIZE/2+size))</x>#g" ./xScale/config_"${type}"_XS"${size}"_1.xml
    sed -i -r "s#<outputprefix>CP_BI_XS_${type}</outputprefix>#<outputprefix>CP_BI_XS${size}_${type}_0</outputprefix>#g" ./xScale/config_"${type}"_XS"${size}"_0.xml
    sed -i -r "s#<outputprefix>CP_BI_XS_${type}</outputprefix>#<outputprefix>CP_BI_XS${size}_${type}_1</outputprefix>#g" ./xScale/config_"${type}"_XS"${size}"_1.xml

    rm -f ./xScale/config_"${type}"_XS"${size}"_run.xml
    cp ./xScale/config_"${type}"_run.xml ./xScale/config_"${type}"_XS"${size}"_run.xml
    sed -Ei "s#<lx unit=\"reduced\">.*</lx>#<lx unit=\"reduced\">$((MSIZE+2*size))\.</lx>#g" ./xScale/config_"${type}"_XS"${size}"_run.xml
    sed -Ei "s#<header>CP_BI_XS#<header>CP_BI_XS${size}#g" ./xScale/config_"${type}"_XS"${size}"_run.xml
    sed -Ei "s#<data>CP_BI_XS#<data>CP_BI_XS${size}#g" ./xScale/config_"${type}"_XS"${size}"_run.xml
    sed -Ei "48s#<x>.*</x>#<x>$((MSIZE/2+size))</x>#" ./xScale/config_"${type}"_XS"${size}"_run.xml
    sed -Ei "60s#<x>.*</x>#<x>$((MSIZE/2+size))</x>#" ./xScale/config_"${type}"_XS"${size}"_run.xml
    sed -Ei "61s#<x>.*</x>#<x>$((MSIZE+2*size))</x>#" ./xScale/config_"${type}"_XS"${size}"_run.xml
    sed -Ei "96s#<lowX>.*</lowX>#<lowX>$((MSIZE/2+size-20))</lowX>#" ./xScale/config_"${type}"_XS"${size}"_run.xml
    sed -Ei "97s#<highX>.*</highX>#<highX>$((MSIZE/2+size+20))</highX>#" ./xScale/config_"${type}"_XS"${size}"_run.xml

    rm -f CI_INIT_XS"${size}"_"${type}"_0.cmd
    rm -f CI_INIT_XS"${size}"_"${type}"_1.cmd

    touch CI_INIT_XS"${size}"_"${type}"_0.cmd
    touch CI_INIT_XS"${size}"_"${type}"_1.cmd
    echo "#!/bin/bash
  #SBATCH -J CI_INIT_XS${size}_${type}_0
  #SBATCH -o ../../../log/%x.out
  #SBATCH -D ../../../build/src
  #SBATCH --partition=medium
  #SBATCH --get-user-env
  #SBATCH --mail-type=all
  #SBATCH --mem=800mb
  #SBATCH --mail-user=hocksa@hsu-hh.de
  #SBATCH --export=NONE
  #SBATCH --time=00:10:00
  #SBATCH --nodes=8
  #SBATCH --ntasks-per-node=4
  #SBATCH --cpus-per-task=18
  #SBATCH --error=../../../log/%x.err

  export OMP_NUM_THREADS=18
  ulimit -s 1000000

  module load slurm_setup
  module load mpi/2021.6.0

  mpiexec -n 32 ./MarDyn ../../examples/AdResS/BlockInterface/xScale/config_${type}_XS${size}_0.xml --loop-abort-time=570
  " > CI_INIT_XS"${size}"_"${type}"_0.cmd
    echo "#!/bin/bash
  #SBATCH -J CI_INIT_XS${size}_${type}_1
  #SBATCH -o ../../../log/%x.out
  #SBATCH -D ../../../build/src
  #SBATCH --partition=medium
  #SBATCH --get-user-env
  #SBATCH --mail-type=all
  #SBATCH --mem=800mb
  #SBATCH --mail-user=hocksa@hsu-hh.de
  #SBATCH --export=NONE
  #SBATCH --time=00:10:00
  #SBATCH --nodes=8
  #SBATCH --ntasks-per-node=4
  #SBATCH --cpus-per-task=18
  #SBATCH --error=../../../log/%x.err

  export OMP_NUM_THREADS=18
  ulimit -s 1000000

  module load slurm_setup
  module load mpi/2021.6.0

  mpiexec -n 32 ./MarDyn ../../examples/AdResS/BlockInterface/xScale/config_${type}_XS${size}_1.xml --loop-abort-time=570
  " > CI_INIT_XS"${size}"_"${type}"_1.cmd
  done
done
