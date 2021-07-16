# Injection scenario

In this scenario, a liquid jet is injected into an evacuated chamber. Therefore, a tube consisting of fixed LJ centers, arranged on dense fcc grid points is generated on the left side of a cuboid system. The liquid is pushed through this tube with constant velocity and replenished by a reservoir to attain stationary conditions. Such a reservoir was also used in a scenario for stationary evaporation across a planar vapor-liquid interface, published by Heinen et al., *J. Chem. Phys.* **151**, 044704 (2019).

## File path structure

A tree view of the file path structure can be found here: [path_info.txt](path_info.txt "File path structure")

## Small system

### Preparation

#### Create liquid reservoir configuration

```bash
cd /path-to-mardyn/examples/Injection/liq/sim01/run03
mpirun -np <numprocs> /path-to-mardyn/src/MarDyn config.xml
```

#### Create dummy reservoir configuration (for the tube)

```bash
cd /path-to-mardyn/examples/Injection/liq/sim01/run04
mpirun -np <numprocs> /path-to-mardyn/src/MarDyn config.xml
```

#### Create scenario start configuration

```bash
cd /path-to-mardyn/examples/Injection/nemd/sim01/run01
mpirun -np <numprocs> /path-to-mardyn/src/MarDyn config.xml
```

### Simulation

#### Start NEMD simulation

```bash
cd /path-to-mardyn/examples/Injection/nemd/sim01/run02
mpirun -np <numprocs> /path-to-mardyn/src/MarDyn config.xml
```
