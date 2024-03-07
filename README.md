ls1-MarDyn Overview
===================

ls1-MarDyn is a massively parallel Molecular Dynamics (MD) code for large systems. Its main target is the simulation of thermodynamics and nanofluidics. ls1-MarDyn is designed with a focus on performance and easy extensibility.


Getting Started
===============

Documentation
--------------
The current doxygen documentation can be found [here](https://ls1mardyn.github.io/ls1-mardyn/).

Prerequisites
--------------
### mandatory requirements
* a C++17 compiler (GCC, Clang, Intel, PGI, Cray, NEC SX, IBM XL, ...)
* a working MPI installation compatible with the MPI 3.0 specification or later (Open MPI, MPICH, MVAPICH, Intel MPI, Cray MPI, NEC MPI, IBM Platform MPI, ...)

### optional requirements
* [FFTW3](http://www.fftw.org)
* [VTK](http://www.vtk.org)
* [QuickSched](https://arxiv.org/abs/1601.05384)


Installation
------------

### Installing ls1-MarDyn using cmake
This is the recommended way of building ls1-MarDyn.

#### Quick guide

Run the following commands to build ls1-MarDyn with the Clang compiler. Adjust options (e.g. ENABLE_MPI) according to the individual needs.
```bash
mkdir build
cd build
CC=clang CXX=clang++ ccmake ..
make -j $(nproc)
```

#### Detailed guide
To build ls1-MarDyn using cmake first create an additional directory on the root ls1-MarDyn directory and change into that directory.
```bash
mkdir build
cd build
```
Next, `cmake` has to be executed. In most cases, you will have to specify the compiler with which ls1-MarDyn should be built:
```bash
CC=clang CXX=clang++ cmake ..
```
Some of the most common compilers are
| Compiler      | CC    | CXX     |
|---------------|-------|---------|
| GCC           | `gcc`   | `g++`     |
| Intel oneAPI  | `icx`   | `icpx`    |
| Clang         | `clang` | `clang++` |

The Intel Classic Compiler (`icc` and `icpc`) is not recommended to use, since it is deprecated and does not work with AutoPas.

Specifying the compiler this way is only possible at the first execution of cmake.
If you want to change the compiler later on, either add another build directory, or first clear the existing build directory.

To configure the options within ls1-MarDyn it is recommended to use `ccmake`:
```bash
ccmake ..
```
That way you can easily edit the available options.

Alternatively, specify the configuration with use of the `cmake` command:
```bash
CC=clang CXX=clang++ cmake -DENABLE_MPI=ON ..
```

Finally, build ls1-MarDyn using:
```bash
make
```

For a faster build, you can make use of parallel building:
```bash
make -j $(nproc)
```

The executable is then found at `build/src/MarDyn`.


### Installing ls1-MarDyn using make

ls1-MarDyn is build from source code using GNU make or alternatively using cmake (see below).

A default build using the GNU compiler and a MPI library providing the mpicxx compiler wrapper is done with
```sh
  cd src
  make
```
To get an overview of options to control the build process, e.g. to use another compiler, disable MPI, ... run
```sh
  make help
```
To see a list of all supported target platforms and compilers call
```sh
  make cfg_list
```
and run then make with the desired configuration:
```sh
  make CFG=<config name>
```
To display further information about the available suboptions for a configuration use
```sh
  make CFG=<cfg name> cfg_help
```


#### ADIOS2 support
If a visualisation with MegaMol and ADIOS2 is desired, an installation of ADIOS2 is needed. By default, ADIOS2 is downloaded and built automatically during the build process of ls1-MarDyn. If connections to external resources, e.g. on HPC systems, are blocked, the following steps (for an MPI build) are required to build ls1-MarDyn with ADIOS2.
Download ADIOS2 locally
```bash
git clone https://github.com/ornladios/ADIOS2.git ADIOS2
```
It may be required to checkout the correct version of ADIOS2:
```bash
cd ADIOS2
git checkout v2.7.1.436
```
Upload it together with ls1-MarDyn to the target HPC system. On the HPC system, after loading the proper modules, create a build folder in the ADIOS2 directory and run cmake. Note: The `PATH-TO-ADIOS2` string must be adjusted.
```bash
cd ADIOS2
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=PATH-TO-ADIOS2/ADIOS2/build/install -DADIOS2_BUILD_EXAMPLES=OFF -DADIOS2_BUILD_TESTING=OFF -DADIOS2_USE_Profiling=OFF -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx
```
Build now using:
```bash
make -j install
```
The `-j` parameter can be adjusting to the present system.

After successfully building ADIOS2, ls1-MarDyn can be build.
In accordance to the above installation steps (section "Configuration"), an additional build directory in the root directory of ls1-MarDyn must be created. In this directory run `cmake`. Note: The `PATH-TO-ADIOS2` string must be adjusted as it was done before.
```bash
mkdir build
cd build
cmake .. -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DADIOS2_DIR=PATH-TO-ADIOS2/ADIOS2/build/install/lib64/cmake/adios2 -DENABLE_MPI=ON -DFIND_PACKAGE_ADIOS2=ON
```

Finally, build ls1-MarDyn using:
```bash
make
```
For a parallel and faster build please use `make`'s `-j` parameter with an appropriate number of tasks.

Note: For both, ADIOS2 and ls1-MarDyn, `ccmake` can be used to configure options.

Running ls1-MarDyn
------------------
The basic command to run ls1-MarDyn is as follows:
```sh
MarDyn [options] <inputfile>
```
where `MarDyn` is the executable build in the INSTALLATION section, `[options]` are any "--"-prefixed options as listed by `MarDyn --help` and `<inputfile>` is a input file.

Detailed help can be obtained by running
```sh
MarDyn --help
```
### running examples
ls1-MarDyn comes with a set of examples, which can be found in the examples folder.
```sh
cd examples/EOX/305K_liq
mpirun -np 2 ../../../src/MarDyn config.xml  --steps 10
```
optional: to make the simulation aware of time limits like on a compute node, which should stop the simulation even if the desired amount of steps is not reached use ```loop-abort-time``` in (s) in XML or CLI:
```sh
mpirun -np 2 ../../../src/MarDyn config.xml  --steps 10 --loop-abort-time 3600
```

AutoPas Support
------------------
ls1-MarDyn supports AutoPas as a replacement for the used linked cells container and the built-in force calculation.

### Building for AutoPas 
To enable support for [AutoPas](https://github.com/AutoPas/AutoPas/), you will have to enable the option `ENABLE_AUTOPAS`.

### Running using AutoPas
To use AutoPas a few modifications to the normal `xml` config files have to be performed:
- The `datastructure` section has to be changed to type `AutoPas`.
- If inside of the `datastructure` section no additional information is given AutoPas will run without auto-tuning and a linked cells container (rebuild frequency = 1, skin = 0).
- Multiple further options can be specified for AutoPas.
  For a quick overview check config_autopas_allOptions.xml in the Argon example directory.
  Additional information for the options can be found in the [official documentation](https://www5.in.tum.de/AutoPas/doxygen_doc/master/namespaceautopas_1_1options.html)
  and within the [responsible readXML method](https://www5.in.tum.de/mardyn/doxygen_doc/html/classAutoPasContainer.html).

### Limitations
- Using AutoPas, currently, only single-centered Lennard-Jones interactions are possible.

Visualisation
------------------
The simulations can be visualised by two external programs which requires the inclusion of the corresponding plugin.

### MegaMol
The MegaMol software is developed by VISUS of the University of Stuttgart. Detailed information on how to install it can be found in its [GitHub repo](https://github.com/UniStuttgart-VISUS/megamol). It supports the import of the following two file formats. See `doc/visualisation_MegaMol.dox` for detailed information.

#### mmpld
This is the old file format for visualisation. Use the `MmpldWriter` plugin to write the visualisation files during simulation.

#### ADIOS2
This kind of visualisation files is recommended. Use the `Adios2Writer` plugin to write the visualisation files during simulation. The produced files can also be used for a restart (see `Adios2Reader`).

### Paraview
Read the documentation in `doc/visualisation_paraview.dox` for detailed information.


Additional resources
====================
ls1-MarDyn is documented using doxygen. To build the documentation run
```sh
doxygen Doxyfile
```
It includes information about the following topics
* \ref ls1MarDynInputFiles Mardyn Input Files
* \ref unitTests Unit tests
* \ref visualisation Graphical Simulation Output

as well as the documentation of the source code.

Contact
=======

* <http://www.ls1-mardyn.de> : the official ls1-MarDyn web page.
* <mailto:contact@ls1-mardyn.de> : can be used for general questions around ls1-MarDyn.
* <mailto:ls1-devel@lists.projects.hlrs.de> : can be used to reach the developers of ls1-MarDyn.



