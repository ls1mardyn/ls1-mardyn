## Sigma
- Value from MolMod is:
    *     \sigma / 1 Armstrong
- Our base unit is 1 nm -> need to rescale, s.t. sigma has unit nm
- 1 A = 0.1 nm
- Component FP CO2_XI:
    *     sig = sig_MM = 3.0354 A = 3.0354 * 0.1 nm = 0.30354 nm 
- Component CG CO2_III:
    *     sig = sig_MM = 4.486 A = 4.486 * 0.1 nm = 0.4486 nm

## Epsilon
- Value from MolMod is:
    *     \epsilon / k_b 
- All our calculations are based on real epsilon, not normalized by k_b
- Need to multiply MolMod eps with k_b = R (Gas Constant)
- Component FP CO2_XI:
    * eps = eps_MM * k_b = 125.317 K * 8.3144621e-3 kJ/(mol*K) = 1.041943447 kJ/mol
- component CG CO2_III:
    * eps = eps_MM * k_b = 189.0 K * 8.3144621e-3 kJ/(mol*K) = 1.571433337 kJ/mol

## Model rescaling
- Initial LJ Site-Site distance = 2.1217446 (reduced)
  - in SI: 2.1217446 * sigma = 0.644034356 nm = 6.44034356 A
- Target O-O distance = 2.32 A
- need to scale down by factor 3 approx.
- pos 1.0608723 -> 0.3536241

## Moment of Inertia
- Value from MolMod has unit: Armstrong² * g / mol
- Need to rescale into nm as distance unit
- 1 A = 0.1 nm => 1 A² = 0.01 nm²
- Component FP CO2_XI:
  *     x: I_x = I_MM_x * 0.01 nm² / A² = 49.53105612429 (A²*g/mol) * 0.01 * (nm²/A²) = 0.4953105612429 nm²*g/mol
  *     y: I_y = I_MM_y * 0.01 nm² / A² = 49.53105612429 (A²*g/mol) * 0.01 * (nm²/A²) = 0.4953105612429 nm²*g/mol
  *     z: I_z = I_MM_z * 0.01 nm² / A² = 0 (A²*g/mol) * 0.01 * (nm²/A²) = 0 nm²*g/mol
- Component other = FP CO2_XI

## Quadrupole
- Unit Force: [F] = kJ / (mol * nm)
- Coulomb's law: F(r) = (-1 / (4 * pi * eps0)) * (q_0 * q_1 / r²)
- Base Unit q: [q] = e = 1.602176565e−19 C = 1.602176565e−19 * A * s
- eps0 = 8.8541878128e−12 * (F / m) = 8.8541878128e−12 * (s⁴ * A²) / (kg * m³)
- eps0 = 3.774144967e-1 * e² * mol / (kJ * nm)
- double q2075 = .75 * absqi * absqj;
- F(r) = -dU(r)/dr = (d/dr) * (6 * Q1 * Q2) / (4 * pi * eps0 * r⁵) * orientation(...)
- F(r) = -5 * (6 * Q1 * Q2) / (4 * pi * eps0 * r⁶) * orientation(...)
- In code factor pi and eps0 are not considered
- In code computed: F*(r) = -5 * (6 * Q1 * Q2) / (4 * r*⁶) * ...
- r is inputted reduced -> need it in SI -> rewrite equation in terms of r*
  - F(r) = -5 * (6 * Q1 * Q2) / (4 * pi * eps0 * r*⁶ * sig⁶)
  - F(r) = -5 * (6 * (Q1/sig³) * (Q2/sig³)) / (4 * pi * eps0 * r*⁶)
  - f := 1 / (pi * eps0)
  - f = 1 / (pi * (3.774144967e-1 * e² * mol / (kJ * nm))) = 0.843396025 * (kJ * nm) / (e² * mol)
  - x := sqrt(f)
  - F(r) = -5 * 6 * (x * Q1 / sig³) * (x * Q2 / sig³) / (4 * r*⁶)
  - now is in correct shape, just need to reduce from F to F*
  - y := sqrt(sig / eps) * x = sqrt(sig * f / eps)
  - F*(r) = F(r) * sig / eps = F(r) * sqrt(sig / eps)²
  - F*(r) = -5 * 6 * (y * Q1 / sig³) * (y * Q2 / sig³) / (4 * r*⁶)
  - Q_i is expected to have unit e * nm², but we get D * Armstrong from MolMod
    - D * A = 0.02081943 e * nm * 0.1 nm = 0.002081943 e * nm² 
- y = sqrt(0.30354 * 0.843396025 * nm² / (e² * 1.041943447)) = 0.495680308 * nm/e
- Q_xml = (y * Q_MM / sig³) = 0.495680308 * 0.002081943 * 3.672697533 / 0.30354³ = 0.135521399

## Dielectric Constant
- Assume as 1 for now, since high temp and low density

## Temperature
- Target temp is 1000K
- T* = k_b * T / eps
- T* = 8.3144621e-3 kJ/(mol*K) * 1000K / (1.041943447 kJ/mol) = 8.3144621 / 1.041943447 = 7.97976332  

