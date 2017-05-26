#! /bin/bash
#export VTKINCLUDEPATH=/usr/include/vtk-6.3
#cd tools/gui

tar xfz ScenarioGenerator.tar.gz
Precision="MARDYN_DPDP"

qmake DEFINES+=$Precision DropletGenerator.pro -o Makefile.droplet
qmake DEFINES+=$Precision CubicGridGenerator.pro -o Makefile.cubic
qmake DEFINES+=$Precision AqueousNaClGenerator.pro -o Makefile.aqueous
qmake DEFINES+=$Precision CrystalLatticeGenerator.pro -o Makefile.crystal
qmake DEFINES+=$Precision MS2RSTGenerator.pro -o Makefile.ms2
qmake DEFINES+=$Precision RayleighTaylorGenerator.pro -o Makefile.rayleigh


if [ -e libMardyn* ]; then
rm -r libMardyn*
fi

#wget www5.in.tum.de/mardyn/lastSuccessfulBuild/lib/libMardyn.so.1.0 .
#ln -s libMardyn.so.1.0 libMardyn.so

make -f Makefile.droplet -j2
make -f Makefile.cubic -j2
make -f Makefile.aqueous -j2
make -f Makefile.crystal -j2
make -f Makefile.ms2 -j2
make -f Makefile.rayleigh -j2

