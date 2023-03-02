# Exploding Liquid

## Description
A liquid film (=xz-plane) is placed in an elongated, otherwise empty domain. Due to the high density the liquid instantly "explodes" along the y-axis. As periodic boundaries are used this simulates an infinitely large liquid plane. Because the plane is placed in the middle of the domain the shock-front of the explosion will collide with its counterpart travelling in the opposite direction upon hitting the periodic boundary in y-direction. 

## Steps
Only one simulation. The start is generated from equilibrated checkpoint files vie the `MultiObjectGenerator`:
- input.header.xml
- input.dat

## Variability
- Can be generated with different thicknesses of the initial liquid. This results in different objects forming behind the explosion shock-front ranging from bubbles to elongated schlieren.

## Notes
- Very quickly changing scenario for two reasons:
    1. At the start almost the whole domain is empty. Then due to the explosion everything is filled.
    2. Unpredictable droplet formation, their movement and merging creates very high and shifting load imbalances. This is interesting for MPI-loadbalancing especially along the y-axis.

## Publications
This experiment was used / shown in the following publications:
- Seckler, S., Gratl, F., Heinen, M., Vrabec, J., Bungartz, H. J., & Neumann, P. (2021). AutoPas in ls1 mardyn: Massively parallel particle simulations with node-level auto-tuning. Journal of Computational Science, 50, 101296.
