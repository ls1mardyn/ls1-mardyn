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

## Notes
 - The checkpoint reading process in config_5 is very memory intensive, and a scenario size of 4k x 6k x 4k fails on both HPE HAWK and Supercomputer Fugaku with a single node.
 - When the liquid droplet is first placed into the spherical void, the sphere surface shows some oscillatory behaviour that is not observed in real life experiments. This is because the liquid boundary is not pre-equilibrated with the vapour. A more realistic scenario is being worked on, where a vapour box is created by deleting particles from a liquid box, keeping a liquid droplet in the centre intact.

## Publications
Refer to the following for more information:
 - Vrabec, J., Kedia, G. K., Fuchs, G., & Hasse, H. (2006). Comprehensive study of the vapour–liquid coexistence of the truncated and shifted Lennard–Jones fluid including planar and spherical interface properties. Molecular physics, 104(09), 1509-1527.
 - Seckler, S., Gratl, F., Heinen, M., Vrabec, J., Bungartz, H. J., & Neumann, P. (2021). AutoPas in ls1 mardyn: Massively parallel particle simulations with node-level auto-tuning. Journal of Computational Science, 50, 101296.
 - Heinen, M., Hoffmann, M., Diewald, F., Seckler, S., Langenbach, K., & Vrabec, J. (2022). Droplet coalescence by molecular dynamics and phase-field modeling. Physics of Fluids, 34(4), 042006.