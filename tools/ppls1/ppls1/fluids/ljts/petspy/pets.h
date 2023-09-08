/*===================================================================================================================================
> Fortran implementation of Pertubated Truncated and Shifted Model Fluid (PeTS) EOS for the article
> T. Hitz, S. Jons, M. Heinen, J. Vrabec, C.-D. Munz
> Comparison of Macro- and Microscopic Solutions of the Riemann Problem II. Two-Phase Shock Tube
> Journal of Computational Physics (2020)
>
> Contact: 
> T. Hitz:  hitz@iag.uni-stuttgart.de
> S. Joens: steven.joens@iag.uni-stuttgart.de
> Institute of Aerodynamics and Gas Dynamics (IAG)
> University of Stuttgart
> Pfaffenwaldring 21
> 70569 Stuttgart
> Germany
> 
> The PeTS EOS was originally published in
> 
> M. Heier, S. Stephan, J. Liu, W. G. Chapman, H. Hasse, K. Langenbach, 
> Equation of state for the Lennard-Jones truncated and shifted fluid with a cut-off radius of 2.5 sigma based on perturbation 
> theory and its applications to interfacial thermodynamics, 
> Molecular Physics 116 (2018) 2083–2094.
> https://doi.org/10.1080/00268976.2018.1447153
===================================================================================================================================*/
/* thermodynamic definitons (from coolprop) */
#define PP_dUMASS_dP_DMASS                        10        /* mol/m^3 UMASS derivated by P at constant DMASS                    */
#define PP_DMOLAR                                 11        /* mol/m^3 Molar density                                             */
#define PP_DMASS                                  12        /* kg/m^3  Mass density                                              */
#define PP_HMOLAR                                 13        /* J/mol   Molar specific enthalpy                                   */
#define PP_HMASS                                  14        /* J/kg    Mass specific enthalpy                                    */
#define PP_P                                      15        /* Pa      Pressure                                                  */
#define PP_Q                                      16        /* mol/mol Mass vapor quality                                        */
#define PP_SMOLAR                                 17        /* J/mol/K Molar specific entropy                                    */
#define PP_SMASS                                  18        /* J/kg/K  Mass specific entropy                                     */
#define PP_T                                      19        /* K       Temperature                                               */
#define PP_UMOLAR                                 20        /* J/mol   Molar specific internal energy                            */
#define PP_UMASS                                  21        /* J/kg    Mass specific internal energy                             */
#define PP_ACENTRIC                               22        /*         Acentric factor                                           */
#define PP_ALPHA0                                 23        /*         Ideal Helmholtz energy                                    */
#define PP_ALPHAR                                 24        /*         Residual Helmholtz energy                                 */
#define PP_A                                      25        /* m/s     Speed of sound                                            */
#define PP_BVIRIAL                                26        /*         Second virial coefficient                                 */
#define PP_CONDUCTIVITY                           27        /* W/m/K   Thermal conductivity                                      */
#define PP_CP0MASS                                28        /* J/kg/K  Ideal gas mass specific constant pressure specific heat   */
#define PP_CP0MOLAR                               29        /* J/mol/K Ideal gas molar specific constant pressure specific heat  */
#define PP_CPMOLAR                                30        /* J/mol/K Molar specific constant pressure specific heat            */
#define PP_CVIRIAL                                31        /*         Third virial coefficient                                  */
#define PP_CVMASS                                 32        /* J/kg/K  Mass specific constant volume specific heat               */
#define PP_CVMOLAR                                33        /* J/mol/K Molar specific constant volume specific heat              */
#define PP_CPMASS                                 34        /* J/kg/K  Mass specific constant pressure specific heat             */
#define PP_DALPHA0_DDELTA_CONSTTAU                35        /*         Derivative of ideal Helmholtz energy with delta           */
#define PP_DALPHA0_DTAU_CONSTDELTA                36        /*         Derivative of ideal Helmholtz energy with tau             */
#define PP_DALPHAR_DDELTA_CONSTTAU                37        /*         Derivative of residual Helmholtz energy with delta        */
#define PP_DALPHAR_DTAU_CONSTDELTA                38        /*         Derivative of residual Helmholtz energy with tau          */
#define PP_DBVIRIAL_DT                            39        /*         Derivative of second virial coefficient with respect to T */
#define PP_DCVIRIAL_DT                            40        /*         Derivative of third virial coefficient with respect to T  */
#define PP_DIPOLE_MOMENT                          41        /*         Dipole moment                                             */
#define PP_FH                                     42        /*         Flammability hazard                                       */
#define PP_FRACTION_MAX                           43        /* mole    maximum value for incompressible solutions                */
#define PP_FRACTION_MIN                           44        /* mole    minimum value for incompressible solutions                */
#define PP_FUNDAMENTAL_DERIVATIVE_OF_GAS_DYNAMICS 45        /*         Fundamental derivative of gas dynamics                    */
#define PP_GAS_CONSTANT                           46        /* J/mol/K Molar gas constant                                        */
#define PP_GMOLAR                                 47        /* J/mol   Molar specific Gibbs energy                               */
#define PP_GWP100                                 48        /*         100-year global warming potential                         */
#define PP_GWP20                                  49        /*         20-year global warming potential                          */
#define PP_GWP500                                 50        /*         500-year global warming potential                         */
#define PP_GMASS                                  51        /* J/kg    Mass specific Gibbs energy                                */
#define PP_GMASS_MS2                              512       /* J/kg    Chemical potential in accordance with ms2                     */
#define PP_HELMHOLTZMASS                          52        /* J/kg    Mass specific Helmholtz energy                            */
#define PP_HELMHOLTZMOLAR                         53        /* J/mol   Molar specific Helmholtz energy                           */
#define PP_HH                                     54        /*         Health hazard                                             */
#define PP_ISOBARIC_EXPANSION_COEFFICIENT         55        /* 1/K     Isobaric expansion coefficient                            */
#define PP_ISOTHERMAL_COMPRESSIBILITY             56        /* 1/Pa    Isothermal compressibility                                */
#define PP_I                                      57        /* N/m     Surface tension                                           */
#define PP_MOLARMASS                              58        /* kg/mol  Molar mass                                                */
#define PP_ODP                                    59        /*         Ozone depletion potential                                 */
#define PP_PCRIT                                  60        /* Pa      Pressure at the critical point                            */
#define PP_PHASE                                  61        /*         Phase index as a float                                    */
#define PP_PH                                     62        /*         Physical hazard                                           */
#define PP_PIP                                    63        /*         Phase identification parameter                            */
#define PP_PMAX                                   64        /* Pa      Maximum pressure limit                                    */
#define PP_PMIN                                   65        /* Pa      Minimum pressure limit                                    */
#define PP_PRANDTL                                66        /*         Prandtl number                                            */
#define PP_PTRIPLE                                67        /* Pa      Pressure at the triple point (pure only)                  */
#define PP_P_REDUCING                             68        /* Pa      Pressure at the reducing point                            */
#define PP_RHOCRIT                                69        /* kg/m^3  Mass density at critical point                            */
#define PP_RHOMASS_REDUCING                       70        /* kg/m^3  Mass density at reducing point                            */
#define PP_RHOMOLAR_CRITICAL                      71        /* mol/m^3 Molar density at critical point                           */
#define PP_RHOMOLAR_REDUCING                      72        /* mol/m^3 Molar density at reducing point                           */
#define PP_SMOLAR_RESIDUAL                        73        /* J/mol/K Residual molar entropy (sr/R = tau*dar_dtau-ar)           */
#define PP_TCRIT                                  74        /* K       Temperature at the critical point                         */
#define PP_TMAX                                   75        /* K       Maximum temperature limit                                 */
#define PP_TMIN                                   76        /* K       Minimum temperature limit                                 */
#define PP_TTRIPLE                                77        /* K       Temperature at the triple point                           */
#define PP_T_FREEZE                               78        /* K       Freezing temperature for incompressible solutions         */
#define PP_T_REDUCING                             79        /* K       Temperature at the reducing point                         */
#define PP_VISCOSITY                              80        /* Pa s    Viscosity                                                 */
#define PP_Z                                      81        /*         Compressibility factor                                    */
#define PP_DIFFUSIVITY                            82        /*                                                                   */
#define PP_DCRIT                                  83        /*                                                                   */
#define PP_GAMMA                                  84        /*                                                                   */
#define PP_QHEAT                                  85        /*                                                                   */
#define PP_PINF                                   86        /*                                                                   */
#define PP_PSAT                                   87        /*                                                                   */
#define PP_H0                                     88        /* J/kg    Enthalpy of formation                                     */
#define PP_S0                                     89        /* J/kg/K  Entropy of formation                                      */
#define PP_dP_dUVOL                               90        /*         P derivated by UVOL                                       */ 
#define PP_dP_dUMASS_DMASS                        91        /*         P derivated by UMASS at constant DMASS                    */
#define PP_dP_dDMASS_UMASS                        92        /*         P derivated by DMASS at constant UMASS                    */
#define PP_RHO_SPNDL                              93        /* kg/m^3  density of the spinodal                                   */
#define PP_A_WOOD                                 94        /* m/s     Speed of sound with Wood approx. in the 2 phase region    */
#define PP_A_THORADE                              95        /* m/s     Speed of sound with Thorade approx. in the 2 phase region */
#define PP_TEST1                                  96        /*         Generic test 1                                            */
#define PP_TEST2                                  97        /*         Generic test 2                                            */
#define PP_TEST3                                  98        /*         Generic test 3                                            */
#define PP_TEST4                                  99        /*         Generic test 4                                            */
#define PP_TEST5                                  100       /*         Generic test 5                                            */
#define PP_TEST6                                  101       /*         Generic test 6                                            */
#define PP_TEST7                                  102       /*         Generic test 7                                            */
#define PP_TEST8                                  103       /*         Generic test 8                                            */
#define PP_TEST9                                  104       /*         Generic test 9                                            */
#define PP_TEST10                                 105       /*         Generic test 10                                           */
#define PP_CONVEX                                 106       /*         Flag indicating convex and concave regions in eos         */
#define PP_dT_dUMASS_DMASS                        107       /*         T derivated by UMASS at constant DMASS                    */
#define PP_dT_dDMASS_UMASS                        108       /*         T derivated by DMASS at constant UMASS                    */
#define PP_CPMASS_HEM                             109       /*         HEM Thorade approximation of cp                           */
#define PP_CVMASS_HEM                             110       /*         HEM Thorade approximation of cv                           */
#define PP_dP_dDMASS_UMASS_HEM                    111       /*         HEM Thorade approximation of Gradient                     */
#define PP_dP_dUMASS_DMASS_HEM                    112       /*         HEM Thorade approximation of Gradient                     */
#define PP_dUMASS_dP_DMASS_HEM                    113       /*         HEM Thorade approximation of Gradient                     */
#define PP_CPMASS_MIX                             114       /*         HEM Mixture approximation of cp                           */
#define PP_CVMASS_MIX                             115       /*         HEM Mixture approximation of cv                           */
#define PP_CONDUCTIVITY_MIX                       116       /*         HEM Mixture approximation of conductiovity                */
#define PP_VISCOSITY_MIX                          117       /*         HEM Mixture approximation of viscosity                    */
#define PP_dP_dDMASS_UMASS_MIX                    118       /*         HEM Mixture approximation of Gradient                     */
#define PP_dP_dUMASS_DMASS_MIX                    119       /*         HEM Mixture approximation of Gradient                     */
#define PP_dUMASS_dP_DMASS_MIX                    120       /*         HEM Mixture approximation of Gradient                     */
#define PP_QTILDEHEAT                             121       /*         Stiff gas parameter                                       */
#define PP_SCRIT                                  122       /*         Critical Value                                            */
#define PP_HCRIT                                  123       /*         Critical Value                                            */
#define PP_ACRIT                                  124       /*         Critical Value                                            */
#define PP_GCRIT                                  125       /*         Critical Value                                            */
#define PP_UCRIT                                  126       /*         Critical Value                                            */
#define PP_DMIN                                   127       /*         min Value                                                 */
#define PP_SMIN                                   128       /*         min Value                                                 */
#define PP_HMIN                                   129       /*         min Value                                                 */
#define PP_AMIN                                   130       /*         min Value                                                 */
#define PP_GMIN                                   131       /*         min Value                                                 */
#define PP_UMIN                                   132       /*         min Value                                                 */
#define PP_QMIN                                   133       /*         min Value                                                 */
#define PP_DMAX                                   134       /*         max Value                                                 */
#define PP_SMAX                                   135       /*         max Value                                                 */
#define PP_HMAX                                   136       /*         max Value                                                 */
#define PP_AMAX                                   137       /*         max Value                                                 */
#define PP_GMAX                                   138       /*         max Value                                                 */
#define PP_UMAX                                   139       /*         max Value                                                 */
#define PP_QMAX                                   140       /*         max Value                                                 */
#define PP_PMID                                   141       /*         mid Value                                                 */
#define PP_TMID                                   142       /*         mid Value                                                 */
#define PP_DMID                                   143       /*         mid Value                                                 */
#define PP_SMID                                   144       /*         mid Value                                                 */
#define PP_HMID                                   145       /*         mid Value                                                 */
#define PP_AMID                                   146       /*         mid Value                                                 */
#define PP_GMID                                   147       /*         mid Value                                                 */
#define PP_UMID                                   148       /*         mid Value                                                 */
#define PP_QMID                                   149       /*         mid Value                                                 */
#define PP_P_SPNDL                                150       /*         spinodal pressure                                         */

/* thermodynamic pair definitons (from coolprop) */
#define PP_DMASS_UMASS                            012021
#define PP_UMASS_DMASS                            021012
#define PP_DMASS_T                                012019
#define PP_T_DMASS                                019012
#define PP_DMASS_P                                012015
#define PP_P_DMASS                                015012
#define PP_DMASS_HMASS                            012014
#define PP_HMASS_DMASS                            014012
#define PP_SMASS_T                                018019
#define PP_T_SMASS                                019018
#define PP_SMASS_P                                018015
#define PP_P_SMASS                                015018
#define PP_SMASS_HMASS                            018014
#define PP_HMASS_SMASS                            014018
#define PP_P_T                                    015019
#define PP_T_P                                    019015
#define PP_Q_T                                    016019
#define PP_T_Q                                    019016
#define PP_Q_P                                    016015
#define PP_P_Q                                    015016
#define PP_A_P                                    025015
#define PP_P_A                                    015025
#define PP_DMASS_SMASS                            012018
#define PP_SMASS_DMASS                            018012
