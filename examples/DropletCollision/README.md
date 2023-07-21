# Droplet Collision

## Description
Two spherical liquid droplets (LJTS) are set in vapor next to each other. A macroscopic velocity was assigned to them so that they will collide eventually.

## Setup
The provided script will setup the initial checkpoint including the two droplets. To make the script run, make sure that the ls1 python toolbox [ppls1](https://github.com/ls1mardyn/ls1-mardyn/tree/master/tools/ppls1/ppls1) works properly. Change the paths in the python script according to the present system.

## Steps
The script will execute three steps, represented in the specified working folder as `s1_Bulk-Liquid`, `s2_Droplet-Equi`, and `s3_Droplet_Collision`.
1. Equilibration of a bulk liquid
2. Bulk liquid is replicated and a droplet is cut out of it. An equilibration simulation for the droplet is executed. Make sure, the steps of equilibration were sufficient (check `Upot` over timesteps).
3. The droplet is duplicated and a macroscopic velocity is assigned respectively.
4. The production simulation as to be started manually from the folder `s3_Droplet_Collision`. Before, the whole simulation folder could be copied to a supercomputer.

## Variability
Just change the parameters like radius or velocity in the python script. The script should take care about the rest.

