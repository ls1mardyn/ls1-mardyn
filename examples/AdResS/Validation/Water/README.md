## Expected Pressure:
- Using van der Waals Equation to estimate pressure 
- We have supercritical, low density water
- Parameters: (https://en.wikipedia.org/wiki/Van_der_Waals_constants_(data_page))
  - a = 0.5536 m⁶*Pa/mol²
  - b = 0.03049 * 1e-3 m³/mol
- Equation: p = RT / (v - b) - a/v²
  - v := V*N_A/N
  - R := Gas Constant
  - T := Temperature
  - V := Volume
  - N_A := Avogardo Constant
  - N := Number Molecules
- Calculation:
  -     V = 36*20*20 nm³ = 14400 nm³ = 14400 * 1e-27 m³ 
  -     N_A = 6.02214076 * 1e+23 / mol
  -     N = 1512
  -     v = V*N_A/N = 14400 * 1e-27 m³ * 6.02214076 * 1e+23 / (mol * 1512) = 8.671882694 m³ / (mol * 1512) = 5.735372152 * 1e-3 * (m³/mol)
  -     p = RT / (v - b) - a/v² = 8.3144621e-3 kJ/(mol*K) * 1000K / (5.735372152 * 1e-3 * (m³/mol) - 0.03049 * 1e-3 m³/mol) - (0.5536 m⁶*Pa/mol²) / (5.735372152 * 1e-3 * (m³/mol))²
  -     p = (8.3144621 * kJ / mol) / (5.704882152 * 1e-3 * (m³/mol)) - (0.5536 / (32.89449372 * 1e-6)) Pa
  -     p = 1.457429247 MPa - 0.01682956439 MPa = 1.440599683 MPa
- Sanity Check using Ideal Gas equation:
  - p* = N * T* / V* = 1512 * 4.805851067 / 14400 = 0.504614362
  - p = p* * eps / sig³ = 0.504614362 * (1.730070696 kJ/mol) / (0.311831716 nm)³ = 28.79138806 * kJ / (mol * nm³) = 28.79138806 * 1.660538921 MPa = 47.80922947 MPa 

## Dielectric Constant
- TODO: need to find more info on this
- Computation based on: 10.1073/pnas.80.14.4575
- T Target is 1000K
- Density: 
  -     Total Volume = 36*20*20 nm³ = 14400 nm³ = 14400 * 1e-27 m³ = 14400 * 1e-21 cm³ = 1.44e-17 cm³  
  -     Total Mass = 1512 * 18.015u = 1512 * 18.015 * 1.660538921 * 1e-27 kg = 4.52308883e-23 kg = 4.52308883e-20 g
  -     Density = Tot_M / Tot_V = 4.52308883e-20 g / 1.44e-17 cm³ = 3.14103391e-3 g/cm³

## Epsilon
- Value from MolMod is:
  *     \epsilon / k_b 
- All our calculations are based on real epsilon, not normalized by k_b
- Need to multiply MolMod eps with k_b = R (Gas Constant)
- Component FP:
  * eps = eps_MM * k_b = 208.07969 K * 8.3144621e-3 kJ/(mol*K) = 1.730070696 kJ/mol
- component CG:
  * eps = eps_MM * k_b = 809.1 K * 8.3144621e-3 kJ/(mol*K) = 6.727231285 kJ/mol

## Sigma
- Value from MolMod is:
  *     \sigma / 1 Armstrong
- Our base unit is 1 nm -> need to rescale, s.t. sigma has unit nm
- 1 A = 0.1 nm
- Component FP:
  *     sig = sig_MM = 3.11831716 A = 3.11831716 * 0.1 nm = 0.311831716 nm 
- Component CG:
  *     sig = sig_MM = 2.641 A = 2.641 * 0.1 nm = 0.2641 nm

## Charges
- Unit Force: [F] = kJ / (mol * nm)
- Coulomb's law: F(r) = (-1 / (4 * pi * eps0)) * (q_0 * q_1 / r²)
- Base Unit q: [q] = e = 1.602176565e−19 C = 1.602176565e−19 * A * s 
- eps0 = 8.8541878128e−12 * (F / m) = 8.8541878128e−12 * (s⁴ * A²) / (kg * m³)
- Sanity Check to see what SI unit we get:
  - F(r) =  (-1 / (4 * pi * 8.8541878128e−12 * (s⁴ * A²) / (kg * m³))) * ((1.602176565e−19)² * A² * s² / nm²)
  - F(r) =  (-(kg * m³) / (1.112650055e-10 * (s⁴ * A²))) * ((1.602176565e−19)² * A² * s² / nm²)
  - F(r) =  -2.307077355e-28 * (kg * m³ * A² * s²) / (s⁴ * A² * nm²)
  - F(r) =  -2.307077355e-31 * (kJ * m) / nm²
  - F(r) =  -2.307077355e-22 * kJ / nm
  - F(r) =  -1.389354458e+2 * kJ / (nm * mol)
- Factor 1/(4 * pi * eps0) is not included in computation in ls1-mardyn
  - Need to split it up and delegate it to the individual charges
- When calculating F(r) r is the only input that is reduced, the rest still has SI units
  - F* = F * sig / eps
  - r* = r / sig
  - F(r) = (1 / (4 * pi * eps0)) * (q_0 * q_1 / r²)    |     express in terms of r* 
  - F(r) = (1 / (4 * pi * eps0)) * (q_0 * q_1 / (r*² * sig²))
  - F(r) = (1 / (4 * pi * eps0)) * ((q_0/sig) * (q_1/sig) / r*²)
  - x := sqrt(1 / (4 * pi * eps0))
  - F(r) = (x * q_0 / sig) * (x * q_1 / sig) / r*²
  - This expression is now in SI units -> need to reduce using F* = F * sig / eps
    - F* = F * sqrt(sig / eps) * sqrt(sig / eps)
    - y := sqrt(sig / eps) * x = sqrt(sig / (4 * pi * eps0 * eps))
    - F(r) = (y * q_0 / sig) * (y * q_1 / sig) / r*²
    - F(r) = (sqrt(sig / (4 * pi * eps0 * eps)) * q_0 / sig) * (sqrt(sig / (4 * pi * eps0 * eps)) * q_1 / sig) / r*²
  - This is now unit-less and the two nominator values in brackets can be used as "charge" values in the xml file
  -     y = sqrt(sig / (4 * pi * eps0 * eps)) = sqrt(0.311831716 / (4 * pi * eps0 * 1.730070696)) = 40248.42779
  -     Charge Positive: y * q+ / sig = 40248.42779 * 0.419547905997027 / 0.311831716 = 54151.46289
  -     Charge Negative: y * q- / sig = 40248.42779 * (-0.839095811994054) / 0.311831716 = 108302.9258


- Second Attempt
  - eps0 = 8.8541878128e−12 * (F / m) = 8.8541878128e−12 * (s⁴ * A²) / (kg * m³)
  - eps0 = 8.8541878128e−12 * (s² * A²) / (J * m) = 8.8541878128e−12 * C² / (J * m)
  - eps0 = 8.8541878128e−12 * (1.602176565e+19)² * e² / (1e-3 * kJ * 1e+9 * nm)
  - eps0 = 2.272843224e+21 * e² / (kJ * nm)
  - eps0 = 2.272843224e+21 * e² * mol / (kJ * nm * mol)
  - eps0 = 2.272843224e+21 * e² * mol / (kJ * nm * 6.02214076 * 1e+23)
  - eps0 = 3.774144967e-1 * e² * mol / (kJ * nm)
  - f := 1 / (4 * pi * eps0)
  - f = 1 / (4 * pi * (3.774144967e-1 * e² * mol / (kJ * nm))) = 0.2108490062 * (kJ * nm) / (e² * mol)
- When calculating F(r) r is the only input that is reduced, the rest still has SI units
  - F* = F * sig / eps
  - r* = r / sig
  - F(r) = (1 / (4 * pi * eps0)) * (q_0 * q_1 / r²)    |     express in terms of r*
  - F(r) = f * (q_0 * q_1 / (r*² * sig²))
  - F(r) = f * ((q_0/sig) * (q_1/sig) / r*²)
  - x := sqrt(f)
  - F(r) = (x * q_0 / sig) * (x * q_1 / sig) / r*²
  - This expression is now in SI units -> need to reduce using F* = F * sig / eps
    - F* = F * sqrt(sig / eps) * sqrt(sig / eps)
    - y := sqrt(sig / eps) * x = sqrt(sig * f / eps) 
    - F(r) = (y * q_0 / sig) * (y * q_1 / sig) / r*²
    - F(r) = (sqrt(sig * f / eps) * q_0 / sig) * (sqrt(sig * f / eps) * q_1 / sig) / r*²
  - This is now unit-less and the two nominator values in brackets can be used as "charge" values in the xml file
  -     y = sqrt(sig * f / eps) = sqrt(0.311831716 * 0.2108490062 / 1.730070696) = 0.194945851
  -     Charge Positive: y * q+ / sig = 0.194945851 * 0.419547905997027 / 0.311831716 = 0.262286097
  -     Charge Negative: y * q- / sig = 0.194945851 * (-0.839095811994054) / 0.311831716 = −0.524572193

## Moment of Inertia
- Value from MolMod has unit: Armstrong² * g / mol
- Need to rescale into nm as distance unit
- 1 A = 0.1 nm => 1 A² = 0.01 nm²
- Component FP:
  *     x: I_x = I_MM_x * 0.01 nm² / A² = 0.89401500129182 (A²*g/mol) * 0.01 * (nm²/A²) = 0.0089401500129182 nm²*g/mol
  *     y: I_y = I_MM_y * 0.01 nm² / A² = 1.6797307537884 (A²*g/mol) * 0.01 * (nm²/A²) = 0.016797307537884 nm²*g/mol
  *     z: I_z = I_MM_z * 0.01 nm² / A² = 2.5737457550802 (A²*g/mol) * 0.01 * (nm²/A²) = 0.025737457550802 nm²*g/mol
- Component CG := FP

## Temperature
- Target temp is 1000K
- T* = k_b * T / eps
- T* = 8.3144621e-3 kJ/(mol*K) * 1000K / (1.730070696 kJ/mol) = 8.3144621 / 1.730070696 = 4.805851067  