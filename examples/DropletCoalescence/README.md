# Droplet Coalescence

## Description
Two spherical liquid droplets (of argon) are suspended in vapour next to each other. Particle attractions cause the droplets to attract each other and slowly coalesce into a single droplet. The two droplets initially form a bridge between them, which eventually grows and subsumes the droplets. 

## Steps
1. `config_1_generateLiq.xml` : Create a liquid argon cube and equilibrate it, saving as a checkpoint.
2. `config_2_replicateLiq.xml` : Take the checkpoint from step 1, replicate it to form a larger liquid body, and equilibrate it to smooth out the edges and create a checkpoint.
3. `config_3_generateVap.xml` : Create a vapour argon cube and equilibrate it, saving as a checkpoint.
4. `config_4_replicateVap.xml` : Take the checkpoint from step 3, replicate it to form a larger vapour body, and equilibrate it to smooth out the edges and create a checkpoint.
5. `config_5_droplet.xml` : Generate the final scenario, by copying the checkpoint from step 4 to create a large vapour box, and then creating two spherical voids inside the box and inserting liquid fron step 2's checkpoint. From here the simulation can be run for as long as required.

## Variability
The droplet size can be any desired size, achieved by editing the values in config_5. However the densities for the fluid and vapour will vary for different droplet sizes, so config_1 and config_3 would need to be rerun with new densities. A python script to calculate densities can be found in tools/scenario-helper/droplet_density_calc.py. The script uses formulae established in reference 1.

These are some common values for quick lookup.
| Diameter | Liquid Density | Vapour Density |
|----------|----------------|----------------|
| 50 nm    | 0.0188022      | 0.000515       |
| 100 nm   | 0.018778       | 0.000506       |
| 200 nm   | 0.0187650      | 0.000501       |

## Notes
 - The checkpoint reading process in config_5 is very memory intensive, since all the particles are read into rank 1 of the simulation and hence may cause the node to run out of memory. A scenario size of 4k x 6k x 4k fails on both HPE HAWK and Supercomputer Fugaku with a single node.
 - When the liquid droplet is first placed into the spherical void, the sphere surface shows some oscillatory behaviour that is not observed in real life experiments. This is because the liquid boundary is not pre-equilibrated with the vapour. A more realistic scenario is being worked on, where a vapour box is created by deleting particles from a liquid box, keeping a liquid droplet in the centre intact.

## Publications
Refer to the following for more information:
 - [1] Vrabec, J., Kedia, G. K., Fuchs, G., & Hasse, H. (2006). Comprehensive study of the vapour–liquid coexistence of the truncated and shifted Lennard–Jones fluid including planar and spherical interface properties. Molecular physics, 104(09), 1509-1527.
 - [2] Seckler, S., Gratl, F., Heinen, M., Vrabec, J., Bungartz, H. J., & Neumann, P. (2021). AutoPas in ls1 mardyn: Massively parallel particle simulations with node-level auto-tuning. Journal of Computational Science, 50, 101296.
 - [3] Heinen, M., Hoffmann, M., Diewald, F., Seckler, S., Langenbach, K., & Vrabec, J. (2022). Droplet coalescence by molecular dynamics and phase-field modeling. Physics of Fluids, 34(4), 042006.