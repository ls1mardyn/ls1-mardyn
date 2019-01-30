ls1-MarDyn Overview
===================

ls1-MarDyn is a massively parallel Molecular Dynamics (MD) code for large systems. Its main target is the simulation of thermodynamics and nanofluidics. ls1-MarDyn is designed with a focus on performance and easy extensibility.


Getting Started
===============

Documentation:
--------------
The current doxygen documentation can be found here:
https://www5.in.tum.de/mardyn/doxygen_doc/html/

Prerequisites:
--------------
### mandatory requirements:
* a C++11 compiler (GCC, Clang, Intel, PGI, Cray, NEC SX, IBM XL, ...)
* a working MPI installation compatible with the MPI 3.0 specification or later (Open MPI, MPICH, MVAPICH, Intel MPI, Cray MPI, NEC MPI, IBM Platform MPI, ...)

### optional requirements:
* FFTW3: <http://www.fftw.org>
* VTK: <http://www.vtk.org>
* QuickSched: <https://arxiv.org/abs/1601.05384>


Installation
------------

ls1-MarDyn is build from source code using GNU make.

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

Running ls1-MarDyn
------------------
The basic command to run ls1-mardyn is as follows:
```sh
MarDyn [options] <inputfile>
```
where MarDyn is the executable build in the INSTALLATION section, `[options]` are any "--"-prefixed options as listed by `MarDyn --help` and `<inputfile>` is a input file.

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



