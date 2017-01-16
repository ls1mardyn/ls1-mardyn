export VTKINCLUDEPATH=/usr/include/vtk-6.2
#cd tools/gui

tar xfz ScenarioGenerator.tar.gz

qmake DropletGenerator.pro -o Makefile.droplet
qmake CubicGridGenerator.pro -o Makefile.cubic
qmake AqueousNaClGenerator.pro -o Makefile.aqueous
qmake CrystalLatticeGenerator.pro -o Makefile.crystal
qmake MS2RSTGenerator.pro -o Makefile.ms2
qmake RayleighTaylorGenerator.pro -o Makefile.rayleigh


if [ -e libMardyn* ]; then
rm -r libMardyn*
fi

cp /home/wwwsccs/html/mardyn/lastSuccessfulBuild/lib/libMardyn.so.1.0 .
ln -s libMardyn.so.1.0 libMardyn.so

make -f Makefile.droplet -j2
make -f Makefile.cubic -j2
make -f Makefile.aqueous -j2
make -f Makefile.crystal -j2
make -f Makefile.ms2 -j2
make -f Makefile.rayleigh -j2
