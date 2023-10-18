# Spinodal Decomposition

## Description
An equilibrated liquid phase is created from a grid generator at, subcritical density rho=0.29, critical temperature T=1.40. In a second step, the temperature is dropped to T=0.70 which triggers the actual spinodal decomposition. The liquid starts to contract and forms 3D patterns which form local inhomogeneities.

## Steps
1. `config_1_genGrid.xml`: Create an equilibrated system at non-critical temperature. Store output in checkpoint.
2. `config_2_spinDec.xml`: Restart the checkpoint with lower temperature which triggers the decomposition.

## Variability
Can be changed arbitrary in size but always needs an very thoroughly equilibrated checkpoint (~500k iterations)

## Publications
This experiment was used / shown in the following publications:
- Gratl, F. A., Seckler, S., Tchipev, N., Bungartz, H. J., & Neumann, P. (2019, May). Autopas: Auto-tuning for particle simulations. In 2019 IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW) (pp. 748-757). IEEE. (not perfectly equilibrated)
- Seckler, S., Gratl, F., Heinen, M., Vrabec, J., Bungartz, H. J., & Neumann, P. (2021). AutoPas in ls1 mardyn: Massively parallel particle simulations with node-level auto-tuning. Journal of Computational Science, 50, 101296.

