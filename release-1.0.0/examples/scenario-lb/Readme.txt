
This scenario should be similar to the ones Martin Buchholz used in his dissertation for loadbalancing.

In the middle we have one drop EOX, surrounded by a gas phase. Concerning the densities, I'm not quite sure, but I guess it doesn't matter as long the ratio of liquid to gas is sufficient (here, it is 100).

16k: "small" example with 16.000 particles and 8 cells along each coordinate dimension. Should be ok to test the load balancing.
90k: 90.000 particles; should be pretty like the original scenario Martin used, as he used 10.000 particles per process. Especially, it should be possible to observe the fact that from 8 to 16 processes there is no speedup without load balancing, as the additional processes only work on cells filled with gas.

Todo: play with the cutoff radius to have more cells?


Observations:

* LB causes some overhead (up to a factor 2 for homogeneous distribution 1CLJ,
  see Diss. Buchholz Abb. 6.6):
** communication is more expensive, as the number of neighbours may be
   much higher
** probably, communication of particles to 26 neighbours cannot be 
   done by 3*2 communications (as in the simple DomainDecomposition)?

* Works the better, the higher load imbalance is, i.e.
** The higher the difference in particle density
** The more expensive the force calculation

* Weak scaling of LB is very good
* Strong scaling is rather poor (on the other hand, load balancing 
  doesn't make sense when there's only one or two cells per process)


Remkarks on implementation:
* as stated in Diss Buchholz, p. 167, number of force calculations 
  could be neglected (alpha is 1.0 anyway, at the moment)
  -> this will reduce the overhead by one traversal of all particle
     pairs
* extra exchange of number of particles and the particles themselves
  -> could be communicated together, as in DomainDecomposition
* don't determine the number of neighbouring processes each time 
  particles are exchanged, but only after rebalancing?
* use SFCDecomposition, could it show better strong scaling?

