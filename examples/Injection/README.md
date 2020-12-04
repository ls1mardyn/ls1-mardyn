# Injection scenario

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
