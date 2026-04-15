# ls1-MarDyn

[![CodeQL](https://github.com/ls1mardyn/ls1-mardyn/actions/workflows/codeql.yml/badge.svg)](https://github.com/ls1mardyn/ls1-mardyn/actions/workflows/codeql.yml)
[![Documentation](https://img.shields.io/badge/docs-doxygen-blue)](https://ls1mardyn.github.io/ls1-mardyn/)
[![License](https://img.shields.io/github/license/ls1mardyn/ls1-mardyn)](https://github.com/ls1mardyn/ls1-mardyn/blob/master/LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Fct500169q-blue)](https://doi.org/10.1021/ct500169q)

ls1-MarDyn is a massively parallel Molecular Dynamics (MD) code designed for large-scale systems. Its primary focus is the simulation of thermodynamics and nanofluidics, with an emphasis on high performance and extensibility.

## Table of Contents

- [Documentation](#documentation)
- [Installation](#installation)
- [Running ls1-MarDyn](#running-ls1-mardyn)
- [AutoPas Support](#autopas-support)
- [Visualisation](#visualisation)
- [Contact](#contact)

## Documentation
The current doxygen documentation can be found [here](https://ls1mardyn.github.io/ls1-mardyn/). It includes information about the following topics
* ls1-MarDyn input files
* Unit tests
* Graphical simulation output
* Source code internals

## Installation

### Prerequisites
#### mandatory requirements
* a C++17 capable compiler
* an MPI installation compatible with the MPI 3.0 specification or later

#### optional requirements
* [Doxygen](https://www.doxygen.nl/) to build the ls1-MarDyn documentation
* [Adios2](https://github.com/ornladios/ADIOS2) for support of reading and writing files in the Adios2 format
* [FFTW3](http://www.fftw.org)
* [Xerces](https://xerces.apache.org/xerces-c/) for support of creating [VTK](http://www.vtk.org) files
* [QuickSched](https://arxiv.org/abs/1601.05384)
* [AutoPas](https://autopas.github.io/doxygen_documentation/git-master/) for alternative particle containers and auto tuning (if enabled, make sure to have the [AutoPas dependencies](https://github.com/AutoPas/AutoPas/blob/master/docs/userdoc/Building.md) installed, as well)

### Installing ls1-MarDyn using cmake
ls1-MarDyn uses CMake as its build system. Run the following commands to configure ls1. Adjust options (e.g. `ENABLE_MPI`) according to your needs.
```sh
git clone https://github.com/ls1mardyn/ls1-mardyn.git
cd ls1-mardyn
cmake -B build/ -S . -DENABLE_MPI=ON
```

ls1-MarDyn provides many configuration options. To list all configuration options run:
```sh
cmake -LAH build/
```
Alternatively use the interactive interface:
```sh
ccmake build/
```

To build the ls1-MarDyn executable run the following command. To speed up the build process it is recommended to run parallel build using the `-j` option:
```sh
cmake --build build/ -j
```
The executable is then found at `build/src/MarDyn`.

To build the documentation of ls1-MarDyn locally using doxygen run
```sh
doxygen Doxyfile
```



## Running ls1-MarDyn
The basic command to run ls1-MarDyn is as follows:
```sh
MarDyn [options] <inputfile>
```
where `MarDyn` is the executable build in the installation section. The available options can be listed with `MarDyn --help`, and `<inputfile>` is an ls1-MarDyn input file.

To understand how to write an input file check out [examples/all-options.xml](https://github.com/ls1mardyn/ls1-mardyn/blob/master/examples/all-options.xml), the various examples in the examples folder and the documentation of the various `readXML()` methods, e.g. via our doxygen documentation.

### Running Examples
ls1-MarDyn comes with a set of examples, which can be found in the examples folder.

```sh
cd examples/EOX/305K_liq
mpirun -np 2 ../../../build/src/MarDyn config.xml  --steps 10
```
optional: to make the simulation aware of time limits like on a compute node, which should stop the simulation even if the desired amount of steps is not reached use ```loop-abort-time``` in (s) in XML or CLI:
```sh
mpirun -np 2 ../../../build/src/MarDyn config.xml  --steps 10 --loop-abort-time 3600
```

## AutoPas Support
ls1-MarDyn supports AutoPas as a replacement for the built-in linked cells container and force calculation.

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

#### Limitations:
- Using AutoPas, currently, only single-centered Lennard-Jones interactions are possible.

## Visualisation
ls1-MarDyn supports visualizing its simulations with various external programs. Specifically it provides plugins for MegaMol and Paraview.

### MegaMol
The MegaMol software is a visualisation software package developed by VISUS of the University of Stuttgart. Detailed information on how to install it can be found in its [GitHub repo](https://github.com/UniStuttgart-VISUS/megamol). Visualisation with MegaMol is typically performed as a post-processing step. Input datasets for MegaMol can be provided via various file formats. See `doc/visualisation_MegaMol.dox` for detailed information.

#### ADIOS2
ADIOS2 is the recommended format to be used with MegaMol. Use the `Adios2Writer` plugin to write the visualisation files during simulation. The produced files can also be used for a restart (see `Adios2Reader`).

#### mmpld
This is the old MegaMol file format. Use the `MmpldWriter` plugin to write the visualisation files during simulation.

### Paraview
Read the documentation in `doc/visualisation_paraview.dox` for detailed information.

# Contact

* [ls1-mardyn.de](http://www.ls1-mardyn.de) - the official ls1-MarDyn web page.
* [general contact](mailto:contact@ls1-mardyn.de) - general questions around ls1-MarDyn
* [developer mailinglist](mailto:ls1-devel@lists.projects.hlrs.de) - mailinglist related to the development of ls1-MarDyn

## Citation
If you use ls1-MarDyn in scientific work, please cite:

```bibtex
@article{niethammer-2014-ls1mardyn,
  author  = {Niethammer, Christoph and Becker, Stefan and Bernreuther, Martin and Buchholz, Martin and Eckhardt, Wolfgang and Heinecke, Alexander and Werth, Stephan and Bungartz, Hans-Joachim and Glass, Colin W. and Hasse, Hans and Vrabec, Jadran and Horsch, Martin},
  title   = {ls1 mardyn: The Massively Parallel Molecular Dynamics Code for Large Systems},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2014},
  volume  = {10},
  number  = {10},
  pages   = {4455--4464},
  doi     = {10.1021/ct500169q}
}
