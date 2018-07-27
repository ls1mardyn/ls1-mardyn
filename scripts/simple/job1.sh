#!/bin/bash
#PBS -N mardyn
#PBS -l nodes=1:ppn=2
#PBS -l walltime=02:00:00
#PBS -o 1nodes.out
#PBS -e 1nodes.err
# Change to the directory that the job was submitted from
cd $PBS_O_WORKDIR
# Launch the parallel job to the allocated compute nodes
#printenv | grep PBS
HYPERTHREADING=1
npes=$PBS_NP
ppn=$PBS_NUM_PPN
omps=$((24*$HYPERTHREADING/$ppn))
nodes=$(($npes/$ppn))
echo ""
echo "running with $ppn processes per node ($omps omp threads per process) on $nodes nodes."
echo "using $HYPERTHREADING threads per core."
echo ""
export OMP_NUM_THREADS=$omps

module list
#aprun -n $npes -N $ppn -d $OMP_NUM_THREADS -cc depth a.out
executable="MarDyn.PAR_RELEASE_AVX2_OPENMP_5044_gcc_NOGEN_noenergy"
aprun -n $npes -N $ppn -d $OMP_NUM_THREADS -cc depth -j $HYPERTHREADING ../$executable 1nodes.xml 50 --final-checkpoint=0
#aprun -n $npes -N $ppn -d $OMP_NUM_THREADS -cc depth -j $HYPERTHREADING ../../exec/$executable ljfluid1node_few.cfg 10 --final-checkpoint=0
#aprun -n 48 -N 24 ./exec/MarDyn.PAR_DEBUG_AOS_4251 inp/simple-lj.cfg 10 --final-checkpoint=0

