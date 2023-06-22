#!/usr/bin/zsh

: "${CLEAN:=0}"
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
    done
  done
  exit 0
fi

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
#SBATCH -o ../log/%x.out
#SBATCH -D ../../..
#SBATCH --partition=medium
#SBATCH --get-user-env
#SBATCH --mail-type=all
#SBATCH --mem=800mb
#SBATCH --mail-user=hocksa@hsu-hh.de
#SBATCH --export=NONE
#SBATCH --time=00:10:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --error=../log/%x.err

export OMP_NUM_THREADS=72
ulimit -s 1000000

module load slurm_setup
module load mpi/2021.6.0

mpiexec -n 8 ./build/src/MarDyn ./examples/AdResS/BlockInterface/${n}k/config_${type}_0.xml --loop-abort-time=570 --final-checkpoint=0
" > CI_INIT_${n}k_"${type}"_0.cmd
  echo "#!/bin/bash
#SBATCH -J CI_INIT_${n}k_${type}_1
#SBATCH -o ../log/%x.out
#SBATCH -D ../../..
#SBATCH --partition=medium
#SBATCH --get-user-env
#SBATCH --mail-type=all
#SBATCH --mem=800mb
#SBATCH --mail-user=hocksa@hsu-hh.de
#SBATCH --export=NONE
#SBATCH --time=00:10:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --error=../log/%x.err

export OMP_NUM_THREADS=72
ulimit -s 1000000

module load slurm_setup
module load mpi/2021.6.0

mpiexec -n 8 ./build/src/MarDyn ./examples/AdResS/BlockInterface/${n}k/config_${type}_1.xml --loop-abort-time=570 --final-checkpoint=0
" > CI_INIT_${n}k_"${type}"_1.cmd
  done
done

