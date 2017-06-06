#!/bin/bash
#export VTKINCLUDEPATH=/usr/include/vtk-6.3
#cd tools/gui

tar xfz ScenarioGenerator.tar.gz
#ExtraDefines="MARDYN_DPDP"
ExtraDefines='MARDYN_SPSP MARDYN_WR=1'


qmake DEFINES+="$ExtraDefines" DropletGenerator.pro -o Makefile.droplet
qmake DEFINES+="$ExtraDefines" CubicGridGenerator.pro -o Makefile.cubic
qmake DEFINES+="$ExtraDefines" AqueousNaClGenerator.pro -o Makefile.aqueous
qmake DEFINES+="$ExtraDefines" CrystalLatticeGenerator.pro -o Makefile.crystal
qmake DEFINES+="$ExtraDefines" MS2RSTGenerator.pro -o Makefile.ms2
qmake DEFINES+="$ExtraDefines" RayleighTaylorGenerator.pro -o Makefile.rayleigh


if [ -e libMardyn* ]; then
rm -r libMardyn*
fi

#wget www5.in.tum.de/mardyn/lastSuccessfulBuild/lib/libMardyn.so.1.0 .
#ln -s libMardyn.so.1.0 libMardyn.so

make -f Makefile.droplet clean
make -f Makefile.cubic clean
make -f Makefile.aqueous clean
make -f Makefile.crystal clean
make -f Makefile.ms2 clean
make -f Makefile.rayleigh clean


make -f Makefile.droplet -j2
make -f Makefile.cubic -j2
make -f Makefile.aqueous -j2
make -f Makefile.crystal -j2
make -f Makefile.ms2 -j2
make -f Makefile.rayleigh -j2

