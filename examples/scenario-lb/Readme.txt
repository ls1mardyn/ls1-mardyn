
This scenario should be similar to the ones Martin Buchholz used in his dissertation for loadbalancing.

In the middle we have one drop EOX, surrounded by a gas phase. Concerning the densities, I'm not quite sure, but I guess it doesn't matter as long the ratio of liquid to gas is sufficient (here, it is 100).

16k: "small" example with 16.000 particles and 8 cells along each coordinate dimension. Should be ok to test the load balancing.
90k: 90.000 particles; should be pretty like the original scenario Martin used, as he used 10.000 particles per process. Especially, it should be possible to observe the fact that from 8 to 16 processes there is no speedup without load balancing, as the additional processes only work on cells filled with gas.

Todo: play with the cutoff radius to have more cells?

