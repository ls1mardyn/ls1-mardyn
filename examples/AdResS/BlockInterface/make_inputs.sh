#!/usr/bin/zsh

: "${CLEAN:=1}"
TYPES=('C6H12' 'CH4')

if [[ $CLEAN -eq 1 ]]
then
  for n in 1 4 16 64 256 1024
  do
    for type in "${TYPES[@]}"
    do
      rm -f ./${n}k/config_"${type}"_0.xml
      rm -f ./${n}k/config_"${type}"_1.xml
      rm -f CI_INIT_${n}k_"${type}"_0.cmd
      rm -f CI_INIT_${n}k_"${type}"_1.cmd
      rm -f ./${n}k/CP_BI_${n}k_"${type}"_0-1.restart.*(N)
      rm -f ./${n}k/CP_BI_${n}k_"${type}"_1-1.restart.*(N)
      rm -f CI_RUN_${n}k_"${type}".cmd
      rm -f CI_BENCH_${n}k_"${type}"_N*.cmd(N)
      rm -f ./${n}k/config_"${type}"_H*.xml(N)
      rm -f CI_BENCH_${n}k_"${type}"_H*_N*.cmd(N)
    done
  done
  exit 0
fi

##
## CREATE EQUILIBRATION CONFIG FILES AND SBATCH COMMANDS
##

for n in 1 4 16 64 256 1024
do
  for type in "${TYPES[@]}"
  do
    rm -f ./${n}k/config_"${type}"_0.xml
    rm -f ./${n}k/config_"${type}"_1.xml
    cp ./${n}k/config_"${type}".xml ./${n}k/config_"${type}"_0.xml
    cp ./${n}k/config_"${type}".xml ./${n}k/config_"${type}"_1.xml
    SIZE=$(sed '27q;d' ./${n}k/config_"${type}".xml | sed -r 's#.*<.*>([0-9]*).*</lx>#\1#g' -)
    sed -i -r "s#<lx unit=\"reduced\">.*</lx>#<lx unit=\"reduced\">$((SIZE/2))\.</lx>#g" ./${n}k/config_"${type}"_0.xml
    sed -i -r "s#<lx unit=\"reduced\">.*</lx>#<lx unit=\"reduced\">$((SIZE/2))\.</lx>#g" ./${n}k/config_"${type}"_1.xml
    sed -i -n '60,81!p' ./${n}k/config_"${type}"_0.xml
    sed -i -n '38,59!p' ./${n}k/config_"${type}"_1.xml
    sed -i -n '80,84!p' ./${n}k/config_"${type}"_0.xml
    sed -i -n '80,84!p' ./${n}k/config_"${type}"_1.xml
    sed -i -n '87,96!p' ./${n}k/config_"${type}"_0.xml
    sed -i -n '87,96!p' ./${n}k/config_"${type}"_1.xml
    sed -i -r "s#<lower> <x>.*</x>#<lower> <x>0</x>#g" ./${n}k/config_"${type}"_1.xml
    sed -i -r "s#<upper> <x>.*</x>#<upper> <x>$((SIZE/2))</x>#g" ./${n}k/config_"${type}"_1.xml
    sed -i -r "s#<outputprefix>CP_BI_${n}k_${type}</outputprefix>#<outputprefix>CP_BI_${n}k_${type}_0</outputprefix>#g" ./${n}k/config_"${type}"_0.xml
    sed -i -r "s#<outputprefix>CP_BI_${n}k_${type}</outputprefix>#<outputprefix>CP_BI_${n}k_${type}_1</outputprefix>#g" ./${n}k/config_"${type}"_1.xml

    rm -f CI_INIT_${n}k_"${type}"_0.cmd
    rm -f CI_INIT_${n}k_"${type}"_1.cmd

    touch CI_INIT_${n}k_"${type}"_0.cmd
    touch CI_INIT_${n}k_"${type}"_1.cmd
    echo "#!/bin/bash
#SBATCH -J CI_INIT_${n}k_${type}_0
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

mpiexec -n 32 ./MarDyn ../../examples/AdResS/BlockInterface/${n}k/config_${type}_0.xml --loop-abort-time=570
" > CI_INIT_${n}k_"${type}"_0.cmd
  echo "#!/bin/bash
#SBATCH -J CI_INIT_${n}k_${type}_1
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

mpiexec -n 32 ./MarDyn ../../examples/AdResS/BlockInterface/${n}k/config_${type}_1.xml --loop-abort-time=570
" > CI_INIT_${n}k_"${type}"_1.cmd
  done
done


for n in 2 3
do
  for type in "${TYPES[@]}"
  do
    rm -f ./1024k${n}/config_"${type}"_0.xml
    rm -f ./1024k${n}/config_"${type}"_1.xml
    cp ./1024k${n}/config_"${type}".xml ./1024k${n}/config_"${type}"_0.xml
    cp ./1024k${n}/config_"${type}".xml ./1024k${n}/config_"${type}"_1.xml
    SIZE=$(sed '27q;d' ./1024k${n}/config_"${type}".xml | sed -r 's#.*<.*>([0-9]*).*</lx>#\1#g' -)
    sed -i -r "s#<lx unit=\"reduced\">.*</lx>#<lx unit=\"reduced\">$((SIZE/2))\.</lx>#g" ./1024k${n}/config_"${type}"_0.xml
    sed -i -r "s#<lx unit=\"reduced\">.*</lx>#<lx unit=\"reduced\">$((SIZE/2))\.</lx>#g" ./1024k${n}/config_"${type}"_1.xml
    sed -i -n '60,81!p' ./1024k${n}/config_"${type}"_0.xml
    sed -i -n '38,59!p' ./1024k${n}/config_"${type}"_1.xml
    sed -i -n '80,84!p' ./1024k${n}/config_"${type}"_0.xml
    sed -i -n '80,84!p' ./1024k${n}/config_"${type}"_1.xml
    sed -i -n '87,96!p' ./1024k${n}/config_"${type}"_0.xml
    sed -i -n '87,96!p' ./1024k${n}/config_"${type}"_1.xml
    sed -i -r "s#<lower> <x>.*</x>#<lower> <x>0</x>#g" ./1024k${n}/config_"${type}"_1.xml
    sed -i -r "s#<upper> <x>.*</x>#<upper> <x>$((SIZE/2))</x>#g" ./1024k${n}/config_"${type}"_1.xml
    sed -i -r "s#<outputprefix>CP_BI_1024k${n}_${type}</outputprefix>#<outputprefix>CP_BI_1024k${n}_${type}_0</outputprefix>#g" ./1024k${n}/config_"${type}"_0.xml
    sed -i -r "s#<outputprefix>CP_BI_1024k${n}_${type}</outputprefix>#<outputprefix>CP_BI_1024k${n}_${type}_1</outputprefix>#g" ./1024k${n}/config_"${type}"_1.xml

    rm -f CI_INIT_1024k${n}_"${type}"_0.cmd
    rm -f CI_INIT_1024k${n}_"${type}"_1.cmd

    touch CI_INIT_1024k${n}_"${type}"_0.cmd
    touch CI_INIT_1024k${n}_"${type}"_1.cmd
    echo "#!/bin/bash
#SBATCH -J CI_INIT_1024k${n}_${type}_0
#SBATCH -o ../../../log/%x.out
#SBATCH -D ../../../build/src
#SBATCH --partition=medium
#SBATCH --get-user-env
#SBATCH --mail-type=all
#SBATCH --mem=800mb
#SBATCH --mail-user=hocksa@hsu-hh.de
#SBATCH --export=NONE
#SBATCH --time=00:20:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18
#SBATCH --error=../../../log/%x.err

export OMP_NUM_THREADS=18
ulimit -s 1000000

module load slurm_setup
module load mpi/2021.6.0

mpiexec -n 32 ./MarDyn ../../examples/AdResS/BlockInterface/1024k${n}/config_${type}_0.xml --loop-abort-time=1170
" > CI_INIT_1024k${n}_"${type}"_0.cmd
  echo "#!/bin/bash
#SBATCH -J CI_INIT_1024k${n}_${type}_1
#SBATCH -o ../../../log/%x.out
#SBATCH -D ../../../build/src
#SBATCH --partition=medium
#SBATCH --get-user-env
#SBATCH --mail-type=all
#SBATCH --mem=800mb
#SBATCH --mail-user=hocksa@hsu-hh.de
#SBATCH --export=NONE
#SBATCH --time=00:20:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18
#SBATCH --error=../../../log/%x.err

export OMP_NUM_THREADS=18
ulimit -s 1000000

module load slurm_setup
module load mpi/2021.6.0

mpiexec -n 32 ./MarDyn ../../examples/AdResS/BlockInterface/1024k${n}/config_${type}_1.xml --loop-abort-time=1170
" > CI_INIT_1024k${n}_"${type}"_1.cmd
  done
done

