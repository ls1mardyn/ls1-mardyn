#! /bin/bash
set -x
#export VTKINCLUDEPATH=/usr/include/vtk-6.3
#cd tools/gui

tar -xzf ScenarioGenerator.tar.gz
Precision="MARDYN_DPDP"
qmake DEFINES+=$Precision DropletGenerator.pro -o Makefile.droplet
qmake DEFINES+=$Precision CubicGridGenerator.pro -o Makefile.cubic
qmake DEFINES+=$Precision AqueousNaClGenerator.pro -o Makefile.aqueous
qmake DEFINES+=$Precision CrystalLatticeGenerator.pro -o Makefile.crystal
qmake DEFINES+=$Precision MS2RSTGenerator.pro -o Makefile.ms2
qmake DEFINES+=$Precision RayleighTaylorGenerator.pro -o Makefile.rayleigh


if [ -f libMardyn.so ]; then
rm libMardyn.so
fi

cp ../../src/libMardyn.so.1.0 libMardyn.so
#ln -s libMardyn.so.1.0 libMardyn.so

if [ ! -f libMardyn.so ]; then
   echo "libMardyn not found. Compile first."
   exit
fi
make -f Makefile.droplet -j2
make -f Makefile.cubic -j2
make -f Makefile.aqueous -j2
make -f Makefile.crystal -j2
make -f Makefile.ms2 -j2
make -f Makefile.rayleigh -j2

