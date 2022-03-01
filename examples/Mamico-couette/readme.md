The two configuration files in this folder are used for preliminary couette tests during the MaMiCo coupling project. However they may be used for testing whenever a generic densely-packed particle box is needed. The example uses the Autopas library and argon molecules.

## Units

The examples use reference units for argon molecules as laid out in the book Numerical Simulation in Molecular Dynamics by Griebel, Knapek and Zumbusch (table 3.9, page 96). The timestep is chosen to be consistent with the timestep used in the coupled MaMiCo simulation, and hence, in the absence of a coupled simulation, may be edited as required.

## Simulation steps

In a coupled simulation, ls1 does not terminate the MD simulation on its own, hence the number of timesteps can be any value.

## Ensemble

Both config files specify a 30 x 30 x 30 box, with 21952 molecules (28 in each spatial direction). The file **ls1configNoCP.xml** builds this box using a cubic generator with density of approximately 0.813. **ls1config.xml** instead uses the file **ls1_CheckpointSimpleMDPeriodic_0__10000_0.inp** as an initial checkpoint. This checkpoint file contains 21952 particles equilibrated for 10000 timesteps using the MD library SimpleMD. Either config file can be used for testing purposes as needed.

## Thermostat

In a coupled simulation, MaMiCo uses its own thermostat on the MD domain, thus the thermostat section in the config files is left empty.

## Parallelization

The config file **ls1config.xml** specifies a domain decomposition of 2 x 2 x 1 for testing. May be edited as needed.

