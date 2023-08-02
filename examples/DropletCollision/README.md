# Droplet Collision

## Description
Two spherical liquid droplets (LJTS) are set in vapor next to each other. A macroscopic velocity was assigned to them so that they will collide eventually.

## Setup
The provided [script](dropColl_MD_setup.py) will setup the initial checkpoint including the two droplets. To make the script run, make sure that the ls1 python toolbox [ppls1](../../tools/ppls1/ppls1) works properly. Change the [paths](dropColl_MD_setup.py#L314-L321) in the python script once according to the present system. In the following, only the settings in the [configuration section](dropColl_MD_setup.py#L292-L295) (e.g. radius) have to be changed. By running the script, all the setup of the main simulation will be done automatically.

## Steps
The script will execute the two initialization steps of the simulation, represented in the specified working folder as `s1_Bulk-Liquid` and `s2_Droplet-Equi`. It will also prepare the main simulation `s3_Droplet_Collision` for execution.
1. Equilibration of a bulk liquid
2. Bulk liquid is replicated and a droplet is cut out of it. An equilibration simulation for the droplet is executed. Make sure, the number of steps of equilibration was sufficient by checking if `Upot` is constant (besides normal fluctuations) over timesteps. Depending on the droplet size, this step can require a lot of RAM. A future workaround could be to use the object generators of ls1 mardyn for setting up the drop.
3. The droplet is duplicated and a macroscopic velocity is assigned respectively.
4. The production simulation has to be started manually from the folder `s3_Droplet_Collision`. Before, the whole simulation folder could be copied to a supercomputer.

## Variability
Just change the parameters like radius or velocity in the python script. The script should take care about the rest.

