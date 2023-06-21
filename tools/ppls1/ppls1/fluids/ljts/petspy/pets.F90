!===================================================================================================================================
!> Fortran implementation of Perturbed Truncated and Shifted Model Fluid (PeTS) EOS for the article
!> T. Hitz, S. Joens, M. Heinen, J. Vrabec, C.-D. Munz
!> Comparison of Macro- and Microscopic Solutions of the Riemann Problem II. Two-Phase Shock Tube
!> Journal of Computational Physics (2020)
!>
!> Contact: 
!> T. Hitz    : hitz@iag.uni-stuttgart.de
!> C.-D. Munz : munz@iag.uni-stuttgart.de
!> Institute of Aerodynamics and Gas Dynamics (IAG)
!> University of Stuttgart
!> Pfaffenwaldring 21
!> 70569 Stuttgart
!> Germany
!> 
!> The PeTS EOS was originally published in
!> 
!> M. Heier, S. Stephan, J. Liu, W. G. Chapman, H. Hasse, K. Langenbach, 
!> Equation of state for the Lennard-Jones truncated and shifted fluid with a cut-off radius of 2.5 sigma based on perturbation 
!> theory and its applications to interfacial thermodynamics, 
!> Molecular Physics 116 (2018) 2083–2094.
!> https://doi.org/10.1080/00268976.2018.1447153
!===================================================================================================================================
#include "pets.h"
MODULE PETS

! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,PARAMETER :: EOSIterationMAX=100
REAL,PARAMETER    :: EOSIterationTOL=1.E-10
REAL,PARAMETER    :: EOSIterationEPS=1.E-12
REAL,PARAMETER    :: pi=ACOS(-1.)
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
REAL                :: Tmin,pmin,rhomin
REAL                :: Tmax,pmax,rhomax
REAL                :: Tcrit,pcrit,rhocrit
REAL                :: R_mass
REAL,DIMENSION(0:6) :: a,b
REAL,DIMENSION(1:2) :: c
REAL                :: ig1,ig2
REAL                :: N_pres(5)         ! coefficients for ancillary equation for saturation pressure
REAL                :: t_a(5)            ! exponents for ancillary equation for saturation pressure
REAL                :: N_liq(4),N_vap(4) ! coefficients for ancillary equation for saturation densities
REAL                :: t_liq(4),t_vap(4) ! exponents for ancillary equation for saturation densities
! coefficients for cp polynomials
REAL                :: T0,rho0

ABSTRACT INTERFACE
  SUBROUTINE EOSIterationFunction(var_iter,var_x,var_y)
  REAL,INTENT(OUT)                 :: var_iter
  REAL,INTENT(IN)                  :: var_x
  REAL,INTENT(IN)                  :: var_y
  END SUBROUTINE EOSIterationFunction
END INTERFACE

ABSTRACT INTERFACE
  SUBROUTINE EOSIterationInitialGuessFunction(guess_min,guess_max,guess,var_x,var_y)
  REAL,INTENT(OUT)                 :: guess_min
  REAL,INTENT(OUT)                 :: guess_max
  REAL,INTENT(OUT)                 :: guess
  REAL,INTENT(IN)                  :: var_x
  REAL,INTENT(IN)                  :: var_y
  END SUBROUTINE EOSIterationInitialGuessFunction
END INTERFACE

PROCEDURE(EOSIterationFunction)            ,POINTER :: EOSIterationTarget
PROCEDURE(EOSIterationFunction)            ,POINTER :: EOSIterationDerivative
PROCEDURE(EOSIterationInitialGuessFunction),POINTER :: EOSIterationInitialGuess

CONTAINS

SUBROUTINE Init()
!==================================================================================================================================
! Initializes the PeTS routines.
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! Argument list declaration
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION 
REAL             :: delta0,tau0
!==================================================================================================================================
! constants
R_mass = 1.

! LJTS critical point
rhocrit   = 0.319   ! only an estimate from LJTS
pcrit     = 6.8/70. ! only a rough estimate
Tcrit     = 1.086   ! only an estimate from LJTS
! ancillary equations
! saturation pressure
N_pres(1) = -6.21
N_pres(2) =  1.5
N_pres(3) = -1.92
N_pres(4) =  2.2
N_pres(5) = -4.76
t_a(1) = 1.00
t_a(2) = 1.50
t_a(3) = 3.25
t_a(4) = 4.85
t_a(5) = 6.63
! saturated liquid density
N_liq(1) =  1.45
N_liq(2) = -0.172
N_liq(3) = -0.298
N_liq(4) =  0.295
t_liq(1) = 0.334
t_liq(2) = 0.667
t_liq(3) = 1.25
t_liq(4) = 1.92
! saturated vapour density
N_vap(1) =  1.59809
N_vap(2) = -0.09975
N_vap(3) = -0.4774
N_vap(4) = -2.33736
t_vap(1) = 1.0
t_vap(2) = 1.5
t_vap(3) = 5.94
t_vap(4) = 0.41452

! Parameters for the PeTS correlation
! a
a(0) = 0.690603404
a(1) = 1.189317012
a(2) = 1.265604153
a(3) =-24.34554201
a(4) = 93.67300357
a(5) =-157.8773415
a(6) = 96.93736697
! b
b(0) = 0.664852128
b(1) = 2.107330790
b(2) =-9.597951213
b(3) =-17.37871193
b(4) = 30.17506222
b(5) = 209.3942909
b(6) =-353.2743581
! c
c(1) = 0.127112544
c(2) = 3.052785558

! reference state for ideal gas part
delta0 = (0.001/0.8)/rhocrit
tau0   = Tcrit/0.8
ig1 = -2.5/tau0
ig2 = 1.5 - LOG(delta0) - 1.5*LOG(tau0)

! calculate critical pressure, based on LJTS estimates for the critical conditions
CALL dT2p_PeTS(pcrit,rhocrit,Tcrit)

! EOS limits
Tmin = 0.6*Tcrit
Tmax = 10.*Tcrit
pmin = 0.
pmax = 70.*pcrit

! calculate density limit
! first, use auxiliary values so that mixture routine can be called
rhomin = 1.E-8
! calculate max density
! NOTE: do not calculate min density, as theoretically 0. is allowed. Use a "practical" value instead
rhomax = 1.67 ! starting value found by try and error
CALL IterateRho(rhomax,pmax,Tmin)


END SUBROUTINE Init


SUBROUTINE Evaluate(n_z,            &
                    integer_x,var_x,       & 
                    integer_y,var_y,       & 
                    integer_z,var_z)
!==================================================================================================================================
! Evaluate PeTS EOS
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! Argument list declaration
! INPUT VARIABLES
INTEGER,INTENT(IN) :: n_z
REAL,INTENT(IN)    :: var_x
REAL,INTENT(IN)    :: var_y
INTEGER,INTENT(IN) :: integer_x
INTEGER,INTENT(IN) :: integer_y 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: var_z(n_z)
INTEGER,INTENT(IN) :: integer_z(n_z) 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION 
INTEGER            :: i_z
REAL               :: var_x_calc
REAL               :: var_y_calc
REAL               :: T,rho
REAL               :: pSat,rho_l,rho_v
!==================================================================================================================================

!----------------------------------------------------------------------------------------------------------------------------------
! No Iteration needed
!----------------------------------------------------------------------------------------------------------------------------------

! EOS=f(DMASS,T)
IF((integer_x.EQ.PP_DMASS.AND.integer_y.EQ.PP_T).OR.(integer_x.EQ.PP_T.AND.integer_y.EQ.PP_DMASS)) THEN
  IF(integer_y.EQ.PP_DMASS) THEN
    T  =var_x
    rho=var_y
  ELSE
    rho=var_x
    T  =var_y
  END IF
  
  DO i_z=1,n_z
    IF(integer_z(i_z).EQ.PP_P) THEN
      ! pressure
      CALL dT2p_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_UMASS) THEN
      ! internal energy
      CALL dT2u_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_A) THEN
      !speed of sound
      CALL dT2a_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_CPMASS) THEN
      ! cp
      CALL dT2cp_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_CVMASS) THEN
      ! cv
      CALL dT2cv_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_SMASS) THEN
      ! entropy
      CALL dT2s_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_HMASS) THEN
      ! enthalpy
      CALL dT2h_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_GMASS) THEN
      ! gibbs energy
      CALL dT2g_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_PIP) THEN
      ! phase indication parameter
      CALL dT2pip_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_TCRIT) THEN
      var_z(i_z) = Tcrit
    ELSEIF(integer_z(i_z).EQ.PP_PCRIT) THEN
      var_z(i_z) = pcrit
    ELSEIF(integer_z(i_z).EQ.PP_Q) THEN
      var_z(i_z) = 0.

    ELSEIF(integer_z(i_z) .EQ. PP_dP_dUMASS_DMASS) THEN
      CALL dpdu_rho_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z) .EQ. PP_dP_dDMASS_UMASS) THEN
      CALL dpdrho_u_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z) .EQ. PP_dT_dUMASS_DMASS) THEN
      CALL dTdrho_u_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z) .EQ. PP_dT_dDMASS_UMASS) THEN
      CALL dTdu_rho_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z) .EQ. PP_I) THEN
          var_z(i_z)=0.
    ELSE
      WRITE(*,*)'Var_x = ',integer_x     
      WRITE(*,*)'Var_y = ',integer_y     
      WRITE(*,*)'Var_z = ',integer_z(i_z)
      WRITE(*,*)'No such output option implemented in PeTS low level interface'
      CALL abort()
    ENDIF
  END DO

!EOS=f(Q,T)
ELSEIF((integer_x.EQ.PP_Q.AND.integer_y.EQ.PP_T).OR.(integer_x.EQ.PP_T.AND.integer_y.EQ.PP_Q)) THEN
  IF(integer_y.EQ.PP_Q) THEN
    var_y_calc=var_x
    var_x_calc=var_y
  ELSE
    var_x_calc=var_x
    var_y_calc=var_y
  END IF
  !IF(integer_z(i_z).EQ.PP_HMASS) THEN
  DO i_z=1,n_z
    IF(integer_z(i_z).EQ.PP_RHO_SPNDL) THEN
      CALL CalcSpinodal_PeTS(rho_l,rho_v,var_y_calc)
      IF (var_x_calc.LT.0.5) THEN
        var_z(i_z) = rho_l
      ELSE
        var_z(i_z) = rho_v
      ENDIF
    ELSEIF(integer_z(i_z).EQ.PP_P) THEN
      CALL IteratePSat(var_y_calc,var_z(i_z),rho_l,rho_v)
    ELSEIF(integer_z(i_z).EQ.PP_DMASS) THEN
      CALL IteratePSat(var_y_calc,pSat,rho_l,rho_v)
      IF (var_x_calc.LT.1.E-8) THEN
        ! saturated liquid
        var_z(i_z) = rho_l
      ELSEIF ((var_x_calc-1.0).LT.1.E-8) THEN
        ! saturated vapour
        var_z(i_z) = rho_v
      ELSE
        var_z(i_z) = -999.E99
        WRITE(*,*)'Dense Gas method not yet implemented in LJTS low level interface',-1,var_x_calc
        CALL abort()
      ENDIF
    ELSEIF(integer_z(i_z).EQ.PP_HMASS) THEN
      CALL IteratePSat(var_y_calc,pSat,rho_l,rho_v)
      IF (var_x_calc.LT.1.E-8) THEN
        ! saturated liquid
        CALL dT2h_PeTS(var_z(i_z),rho_l,var_y_calc)
      ELSEIF ((var_x_calc-1.0).LT.1.E-8) THEN
        ! saturated vapour
        CALL dT2h_PeTS(var_z(i_z),rho_v,var_y_calc)
      ELSE
        var_z(i_z) = -999.E99
        WRITE(*,*)'Dense Gas method not yet implemented in PeTS low level interface',-1,var_x_calc
        CALL abort()
      ENDIF
    ELSE
      WRITE(*,*)'Var_x = ',integer_x     
      WRITE(*,*)'Var_y = ',integer_y     
      WRITE(*,*)'Var_z = ',integer_z(i_z)
      WRITE(*,*)'No such output option implemented in PeTS low level interface'
      CALL Abort()
    ENDIF
  END DO

! EOS=f(Q,P)
ELSEIF((integer_x.EQ.PP_Q.AND.integer_y.EQ.PP_P).OR.(integer_x.EQ.PP_P.AND.integer_y.EQ.PP_Q)) THEN
  IF(integer_y.EQ.PP_Q) THEN
    var_y_calc=var_x
    var_x_calc=var_y
  ELSE
    var_x_calc=var_x
    var_y_calc=var_y
  END IF
  CALL IterateTSat(var_y_calc,T,rho_l,rho_v)
  DO i_z=1,n_z
    IF(integer_z(i_z).EQ.PP_T) THEN
      var_z(i_z) = T
    ELSE
      WRITE(*,*)'Var_x = ',integer_x     
      WRITE(*,*)'Var_y = ',integer_y     
      WRITE(*,*)'Var_z = ',integer_z(i_z)
      WRITE(*,*)'No such output option implemented in PeTS low level interface'
      CALL abort()
    ENDIF
  END DO

  ! EOS=f(p,T)
ELSEIF((integer_x.EQ.PP_P.AND.integer_y.EQ.PP_T).OR.(integer_x.EQ.PP_T.AND.integer_y.EQ.PP_P)) THEN

!----------------------------------------------------------------------------------------------------------------------------------
! Newton Iteration if only p and T are known
!----------------------------------------------------------------------------------------------------------------------------------
  IF(integer_y.EQ.PP_P) THEN
    var_y_calc=var_x
    var_x_calc=var_y
  ELSE
    var_x_calc=var_x
    var_y_calc=var_y
  END IF
  ! iterate rho
  CALL IterateRho(rho,var_x_calc,var_y_calc)
  DO i_z=1,n_z
    IF(integer_z(i_z).EQ.PP_DMASS) THEN
      ! density
      var_z(i_z)=rho
    ELSEIF(integer_z(i_z).EQ.PP_P) THEN
      ! pressure
      CALL dT2p_PeTS(var_z(i_z),rho,var_y_calc)
    ELSEIF(integer_z(i_z).EQ.PP_UMASS) THEN
      ! internal energy
      CALL dT2u_PeTS(var_z(i_z),rho,var_y_calc)
    ELSEIF(integer_z(i_z).EQ.PP_A) THEN
      !speed of sound
      CALL dT2a_PeTS(var_z(i_z),rho,var_y_calc)
    ELSEIF(integer_z(i_z).EQ.PP_CPMASS) THEN
      ! cp
      CALL dT2cp_PeTS(var_z(i_z),rho,var_y_calc)
    ELSEIF(integer_z(i_z).EQ.PP_CVMASS) THEN
      ! cv
      CALL dT2cv_PeTS(var_z(i_z),rho,var_y_calc)
    ELSEIF(integer_z(i_z).EQ.PP_SMASS) THEN
      ! entropy
      CALL dT2s_PeTS(var_z(i_z),rho,var_y_calc)
    ELSEIF(integer_z(i_z).EQ.PP_HMASS) THEN
      ! enthalpy
      CALL dT2h_PeTS(var_z(i_z),rho,var_y_calc)
    ELSEIF(integer_z(i_z).EQ.PP_GMASS) THEN
      ! gibbs energy
      CALL dT2g_PeTS(var_z(i_z),rho,var_y_calc)
    ELSEIF(integer_z(i_z).EQ.PP_PIP) THEN
      ! phase indication parameter
      CALL dT2pip_PeTS(var_z(i_z),rho,var_y_calc)
    ELSEIF(integer_z(i_z).EQ.PP_TCRIT) THEN
      var_z(i_z) = Tcrit
    ELSEIF(integer_z(i_z).EQ.PP_PCRIT) THEN
      var_z(i_z) = pcrit
    ELSEIF(integer_z(i_z) .EQ. PP_dP_dUMASS_DMASS) THEN
      CALL dpdu_rho_PeTS(var_z(i_z),rho,var_y_calc)
    ELSEIF(integer_z(i_z) .EQ. PP_dP_dDMASS_UMASS) THEN
      CALL dpdrho_u_PeTS(var_z(i_z),rho,var_y_calc)
    ELSEIF(integer_z(i_z) .EQ. PP_dT_dUMASS_DMASS) THEN
      CALL dTdrho_u_PeTS(var_z(i_z),rho,var_y_calc)
    ELSEIF(integer_z(i_z) .EQ. PP_dT_dDMASS_UMASS) THEN
      CALL dTdu_rho_PeTS(var_z(i_z),rho,var_y_calc)
    ELSEIF(integer_z(i_z) .EQ. PP_I) THEN
          var_z(i_z)=0.
    ELSE
      WRITE(*,*)'Var_x = ',integer_x     
      WRITE(*,*)'Var_y = ',integer_y     
      WRITE(*,*)'Var_z = ',integer_z(i_z)
      WRITE(*,*)'No such output option implemented in PeTS low level interface'
      CALL abort()
    ENDIF
  END DO

ELSE
!----------------------------------------------------------------------------------------------------------------------------------
! Use iteration routines if density is known
!----------------------------------------------------------------------------------------------------------------------------------

  ! EOS=f(DMASS,U)
  ! PeTS low-level-interface does not support this input directly
  ! *output      variable is always temperature
  ! *iteration   variable is        internal energy 
  ! *calculation variable is        density
  IF((integer_x.EQ.PP_DMASS.AND.integer_y.EQ.PP_UMASS).OR.(integer_x.EQ.PP_UMASS.AND.integer_y.EQ.PP_DMASS)) THEN
    IF(integer_y.EQ.PP_DMASS) THEN
      var_y_calc=var_x
      var_x_calc=var_y
    ELSE
      var_x_calc=var_x
      var_y_calc=var_y
    END IF
    ! choose routines for target function and derivative
    EOSIterationTarget => dT2u_PeTS
    EOSIterationDerivative => dudT_rho_PeTS
    EOSIterationInitialGuess => InitialGuess_InternalEnergy_PeTS
  ! EOS=f(DMASS,P)
  ! PeTS low-level-interface does not support this input directly
  ! *output      variable is always temperature
  ! *iteration   variable is        pressure
  ! *calculation variable is        density
  ELSEIF((integer_x.EQ.PP_DMASS.AND.integer_y.EQ.PP_P).OR.(integer_x.EQ.PP_P.AND.integer_y.EQ.PP_DMASS)) THEN
    IF(integer_y.EQ.PP_DMASS) THEN
      var_y_calc=var_x
      var_x_calc=var_y
    ELSE
      var_x_calc=var_x
      var_y_calc=var_y
    END IF
    ! choose routines for target function and derivative
    EOSIterationTarget => dT2p_PeTS
    EOSIterationDerivative => dpdT_rho_PeTS
    EOSIterationInitialGuess => InitialGuess_VdW_PeTS
  ! EOS=f(DMASS,HMASS)
  ! PeTS low-level-interface does not support this input directly
  ! *output      variable is always temperature
  ! *iteration   variable is        enthalpy
  ! *calculation variable is        density
  ELSEIF((integer_x.EQ.PP_DMASS.AND.integer_y.EQ.PP_HMASS).OR.(integer_x.EQ.PP_HMASS.AND.integer_y.EQ.PP_DMASS)) THEN
    IF(integer_y.EQ.PP_DMASS) THEN
      var_y_calc=var_x
      var_x_calc=var_y
    ELSE
      var_x_calc=var_x
      var_y_calc=var_y
    END IF
    ! choose routines for target function and derivative
    EOSIterationTarget => dT2h_PeTS
    EOSIterationDerivative => dhdT_rho_PeTS
    EOSIterationInitialGuess => InitialGuess_Default_PeTS
  ! EOS=f(DMASS,SMASS)
  ! PeTS low-level-interface does not support this input directly
  ! *output      variable is always temperature
  ! *iteration   variable is        entropy
  ! *calculation variable is        density
  ELSEIF((integer_x.EQ.PP_DMASS.AND.integer_y.EQ.PP_SMASS).OR.(integer_x.EQ.PP_SMASS.AND.integer_y.EQ.PP_DMASS)) THEN
    IF(integer_y.EQ.PP_DMASS) THEN
      var_y_calc=var_x
      var_x_calc=var_y
    ELSE
      var_x_calc=var_x
      var_y_calc=var_y
    END IF
    ! choose routines for target function and derivative
    EOSIterationTarget => dT2s_PeTS
    EOSIterationDerivative => dsdT_rho_PeTS
    EOSIterationInitialGuess => InitialGuess_Entropy_PeTS
  ELSE
    WRITE(*,*)'Var_x = ',integer_x
    WRITE(*,*)'Var_y = ',integer_y
    WRITE(*,*)'No such input pair implemented in PeTS low level interface'
    CALL abort()
  ENDIF

  ! iterate and calculate T
  CALL IterateT(T,var_x_calc,var_y_calc)
  rho = var_x_calc

  DO i_z=1,n_z
    IF(integer_z(i_z).EQ.PP_T) THEN
      ! temperature
      var_z(i_z)=T 
    ELSEIF(integer_z(i_z).EQ.PP_P) THEN
      ! pressure
      CALL dT2p_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_UMASS) THEN
      ! internal energy
      CALL dT2u_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_A) THEN
      !speed of sound
      CALL dT2a_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_CPMASS) THEN
      ! cp
      CALL dT2cp_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_CVMASS) THEN
      ! cv
      CALL dT2cv_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_SMASS) THEN
      ! entropy
      CALL dT2s_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_HMASS) THEN
      ! enthalpy
      CALL dT2h_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_GMASS) THEN
      ! gibbs energy
      CALL dT2g_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z).EQ.PP_PIP) THEN
      ! phase indication parameter
      CALL dT2pip_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z) .EQ. PP_dP_dUMASS_DMASS) THEN
      CALL dpdu_rho_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z) .EQ. PP_dP_dDMASS_UMASS) THEN
      CALL dpdrho_u_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z) .EQ. PP_dT_dUMASS_DMASS) THEN
      CALL dTdrho_u_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z) .EQ. PP_dT_dDMASS_UMASS) THEN
      CALL dTdu_rho_PeTS(var_z(i_z),rho,T)
    ELSEIF(integer_z(i_z) .EQ. PP_I) THEN
        var_z(i_z)=0.
    ELSEIF(integer_z(i_z).EQ.PP_TCRIT) THEN
      var_z(i_z) = Tcrit
    ELSEIF(integer_z(i_z).EQ.PP_PCRIT) THEN
      var_z(i_z) = pcrit
    ELSE
      WRITE(*,*)'Var_x = ',integer_x     
      WRITE(*,*)'Var_y = ',integer_y     
      WRITE(*,*)'Var_z = ',integer_z(i_z)
      WRITE(*,*)'No such output option implemented in PeTS low level interface'
      CALL Abort()
    ENDIF
  END DO

ENDIF


END SUBROUTINE Evaluate


SUBROUTINE IterateT(T,var_x,var_y)
!==================================================================================================================================
!> Iterate the temperature for a given density and a iteration variable, e.g. pressure, internal energy, ...
!> using a Newton algorthim. It uses the exact derivatives of the underlying EOS
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! Argument list declaration
! INPUT VARIABLES
REAL,INTENT(IN) :: var_x
REAL,INTENT(IN) :: var_y
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: T
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION 
REAL    :: var_iter
INTEGER :: i, i_inner
REAL    :: f,dT,res,dvar_iterdT,delta
REAL    :: Tstart
REAL    :: Ta,Tb,Tc
REAL    :: tol
!==================================================================================================================================

! Guess the initial temperature
CALL EOSIterationInitialGuess(Ta,Tb,Tstart,var_x,var_y)

T = Tstart

! Newton = 1 
! set local tolerance
tol = EOSIterationTOL
DO i = 1, EOSIterationMAX
  ! increase the tolerance if Newton takes a while to converge.
  ! this idea is taken from RefProp
  IF (i.EQ.10) tol=10.*tol
  IF (i.EQ.30) tol=10.*tol
  IF (i.EQ.50) tol=10.*tol
  CALL EOSIterationTarget(var_iter,var_x,T)
  f = var_iter-var_y
  res  =  ABS(f)/(ABS(var_y)+1.)
  IF(res.LT.tol) EXIT
  CALL EOSIterationDerivative(dvar_iterdT,var_x,T)
  dT   = -f/dvar_iterdT
  delta= 0.9
  ! Globalization strategy: linesearch   
  DO i_inner=1,10
    Tc = T+delta*dT
    CALL EOSIterationTarget(var_iter,var_x,Tc)
    IF((ABS(var_iter-var_y)/(ABS(var_y)+1.).LT.res)) THEN
      EXIT
    ELSE
      delta = 0.5*delta
    ENDIF
  ENDDO
  T=T+delta*dT
ENDDO

! check if iterationnumber is reached
IF (i.GE.EOSIterationMAX .AND. res.GT.tol) THEN
  WRITE(*,*)'IterateT failed. Tstart',Tstart,'T iterated',T,'for rho',var_x,'and var_y',var_y
  WRITE(*,*)"Temperature could not be iterated in IterateT. i,res",i,ABS(res)
  CALL Abort()
ENDIF

END SUBROUTINE IterateT

SUBROUTINE IterateRho(rho,p,T)
!==================================================================================================================================
!> Iterate the density for a given pressure and temperature and a iteration variable, e.g. pressure, internal energy, ...
!> You can only use the Newton algorithm.
!> The Newton algorithms use the exact derivatives of the underlying EOS
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! Argument list declaration
! INPUT VARIABLES
REAL,INTENT(IN)  :: p
REAL,INTENT(IN)  :: T
REAL,INTENT(OUT) :: rho
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION 
INTEGER          :: i, i_inner
DOUBLE PRECISION :: var_iter
DOUBLE PRECISION :: f,drho,res,dpdrho_T,delta
REAL             :: rho_inner
REAL             :: tol
REAL             :: pSat,rho_l,rho_v
!==================================================================================================================================

! choose an appropriate initial guess for rho
IF(p.GE.pcrit .AND. T.GE.Tcrit) THEN
  ! supercritical fluid
  rho = rhocrit+10. ! don't use the critical density, as it is too close to the two phase regime
ELSEIF(p.LT.pcrit .AND. T.GE.Tcrit) THEN
  ! ideal gas
  rho = p/(R_mass*T)
ELSEIF(T.LT.Tcrit) THEN
  ! subcritical temperature ...
  IF (p.GE.pcrit) THEN
    ! ... but supercritical pressure
    rho = rhomax ! density lies on the liquid side, choose large initial guess
  ELSE
    ! ... but subcritical pressure
    ! the fluid is subcritical in all respects, i.e. there are multiple roots due to the VdW loop
    ! first, calculate the saturation pressure to get a better idea where we are
    CALL IteratePSat(T,pSat,rho_l,rho_v)
    ! in this case, only calculate stable roots, i.e. we return only single phase states
    IF (p.GT.pSat) THEN
      ! pressure is above saturation pressure, we have a stable liquid phase
      rho = rhomax ! choose large initial guess
    ELSEIF(p.LT.pSat) THEN
      ! pressure is below saturation pressure, we have a stable gas phase
      rho = p/(R_mass*T) ! choose ideal gas guess
    ELSE
      ! if the given pressure is the saturation pressure, return liquid density
      rho = rho_l
      RETURN
    ENDIF
  ENDIF
ENDIF

! set local tolerance
tol = EOSIterationTOL
DO i = 1, EOSIterationMAX
  ! increase the tolerance if Netwon takes long to converge.
  ! this idea is taken from RefProp
  IF (i.EQ.10) tol=10.*tol
  IF (i.EQ.30) tol=10.*tol
  IF (i.EQ.50) tol=10.*tol
  CALL dT2p_PeTS(var_iter,rho,T) ! target function
  f = var_iter-p
  res  =  ABS(f)/(ABS(p)+1.)
  IF(res.LT.tol) EXIT
  CALL dpdrho_T_PeTS(dpdrho_T,rho,T) ! derivative of target function
  drho = -f/dpdrho_T
  delta= 0.9
  ! Globalization strategy: linesearch   
  DO i_inner=1,10
    rho_inner = rho+delta*drho
    CALL dT2p_PeTS(var_iter,rho_inner,T) ! target function
    IF((ABS(var_iter-p)/(ABS(p)+1.).LT.res)) THEN
      EXIT
    ELSE
      delta = 0.5*delta
    ENDIF
  ENDDO
  rho=rho+delta*drho
ENDDO

! check if iterationnumber is reached
IF (i.GE.EOSIterationMAX .AND. res.GT.tol) THEN
  WRITE(*,*)"Density could not be iterated in IterateRho. i,res",i,ABS(res)
  CALL Abort()
ENDIF

END SUBROUTINE IterateRho


!==================================================================================================================================
!> Iterates the saturation pressure and densities for a given temperature
!> The procedure follows Akasaka, 2008
!> https://doi.org/10.1299/jtst.3.442
!==================================================================================================================================
SUBROUTINE IteratePSat(T,pSat,rho_l,rho_v)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! Argument list declaration
! INPUT VARIABLES
REAL,INTENT(IN) :: T
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: pSat
REAL,INTENT(OUT) :: rho_l
REAL,INTENT(OUT) :: rho_v
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION 
INTEGER :: i,i_inner
REAL    :: rho_sl,rho_sv
REAL    :: rho(2),rho_ancillary(2)
REAL    :: J(2),K(2)
REAL    :: dJdrho(2),dKdrho(2)
REAL    :: Jdiff,Kdiff
REAL    :: alphar00(2),alphar01(2),alphar02(2)
REAL    :: delta(2),delta_inner(2),tau
REAL    :: ddelta(2)
REAL    :: res(2)
REAL    :: det
REAL    :: gamma
REAL    :: tol
REAL    :: deltaa,deltab,delta_upper(2),deltac,var_m
REAL    :: fa,fb,fc
REAL    :: p,pa,pb,pc
REAL    :: res_p,res_delta_inner,res_g,res_delta
REAL    :: g(2)
!==================================================================================================================================
! NOTE: index 1 - liquid
!       index 2 - vapour

IF (T.GE.Tcrit) THEN
  rho_l = rhocrit
  rho_v = rho_l
  pSat  = pcrit
  RETURN
ENDIF

IF (T.LT.1.E-10) THEN
  rho_l = 0.
  rho_v = rho_l
  pSat  = 0.
  RETURN
ENDIF

!Set initial values
rho(1) = rhomax ! liquid
rho(2) = rhomin ! vapour
CALL AncillaryEquation_DensitySat(T,rho_ancillary)
rho = rho_ancillary

! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit 

tol = 1.E-8
IF(T.LE.0.9*Tcrit)THEN
  !Start Newton-Raphson-Iteration
  DO i = 1, EOSIterationMAX
    IF (i.EQ.10) tol=10.*tol
    IF (i.EQ.30) tol=10.*tol
    IF (i.EQ.50) tol=10.*tol
    ! calculate J and K
    ! J (eq. 21)
    alphar01(1) = ALPHAR(tau,delta(1),0,1)
    alphar01(2) = ALPHAR(tau,delta(2),0,1)
    J = delta*(1.+delta*alphar01)
    ! K (eq. 22)
    alphar00(1) = ALPHAR(tau,delta(1),0,0)
    alphar00(2) = ALPHAR(tau,delta(2),0,0)
    K = delta*alphar01 + alphar00 + LOG(delta)
    ! residual
    Jdiff = J(2)-J(1) ! represents pressure equilibrium
    Kdiff = K(2)-K(1) ! represents Gibbs energy equilibrium

    ! is converged?
    res(1) = Jdiff
    res(2) = Kdiff
    IF (SUM(ABS(res)).LT.tol) EXIT
    
    ! calculate derivative of J (eq.23) and K (eq.24)
    alphar02(1) = ALPHAR(tau,delta(1),0,2)
    alphar02(2) = ALPHAR(tau,delta(2),0,2)
    dJdrho = 1. + 2.*delta*alphar01 + delta*delta*alphar02
    dKdrho = 2.*alphar01 + delta*alphar02 + 1./delta
    ! calculate step size for Newton-Raphson
    det  = dJdrho(2)*dKdrho(1) - dJdrho(1)*dKdrho(2) ! (eq.25)
    ddelta(1) = Kdiff*dJdrho(2) - Jdiff*dKdrho(2) ! (eq.19)
    ddelta(2) = Kdiff*dJdrho(1) - Jdiff*dKdrho(1) ! (eq.20)
    ddelta = ddelta/det

    gamma= 0.9
    ! Globalization strategy: linesearch   
    DO i_inner=1,10
      delta_inner = delta + gamma*ddelta

      ! calculate J and K
      ! J (eq. 21)
      alphar01(1) = ALPHAR(tau,delta_inner(1),0,1)
      alphar01(2) = ALPHAR(tau,delta_inner(2),0,1)
      J = delta_inner*(1.+delta_inner*alphar01)
      ! K (eq. 22)
      alphar00(1) = ALPHAR(tau,delta_inner(1),0,0)
      alphar00(2) = ALPHAR(tau,delta_inner(2),0,0)
      K = delta_inner*alphar01 + alphar00 + LOG(delta_inner)
      ! residual
      res(1) = J(2)-J(1) ! represents pressure equilibrium
      res(2) = K(2)-K(1) ! represents Gibbs energy equilibrium
      !Check the convergence criteria
      IF((ABS(res(1)).LT.ABS(Jdiff)).AND.(ABS(res(2)).LT.ABS(Kdiff)))THEN
        EXIT
      ELSE
        gamma = 0.5*gamma
      ENDIF
    ENDDO
    delta = delta + gamma*ddelta
  END DO
ELSE
  ! Fallback to Bisection in the vicinity of the critical point
  CALL CalcSpinodal_PeTS(rho_sl,rho_sv,T)
  ! Set initial value for the bisection scheme
  !delta(2) = 0.5*(rhomin(1)+rho_sv)/rhocrit
  delta(2) = 0.5*(rho_sv)/rhocrit
  ! Set upper bound of the biseciton scheme
  delta_upper(2) = rho_sv/rhocrit
  ! Starting bisection scheme
  DO i = 1, EOSIterationMAX

    ! Calculate pressure at gaseous side
    CALL dT2p_PeTS(p,delta(2)*rhocrit,T)
    ! Calculation of the density with the same pressure as at the gaseous side
    ! Set initial value for the bisection scheme
    deltaa = rho_sl/rhocrit
    deltab = rhomax/rhocrit
    ! Calculation at the lower interval boundary
    CALL dT2p_PeTS(pa,deltaa*rhocrit,T)
    ! At temperatures near the critical point the pressure of the middle density DMASSc of the bisection scheme
    ! can be lower than the pressure at the liquid spinodale density at the given temperature, therefor a correction
    ! is neccessary, so the density DMASSc has to be increased until the pressure at DMASSc is higher than at liquid
    ! spinodale
    IF(p.LT.pa)THEN
      DO WHILE(p.LT.pa)
        delta_upper(1) = delta(2)
        delta(2) = 0.5*(delta_upper(1)+delta_upper(2))
        CALL dT2p_PeTS(p,delta(2)*rhocrit,T)
      END DO
    END IF
    fa = pa-p
    ! Calculation at the upper interval boundary
    CALL dT2p_PeTS(pb,deltab*rhocrit,T)
    fb = pb-p
    DO i_inner = 1, EOSIterationMAX

      ! Bisection
      deltac = 0.5*(deltaa+deltab) 
      CALL dT2p_PeTS(pc,deltac*rhocrit,T)
      fc = pc-p
      !calculation of residual
      res_p           = fc
      res_delta_inner = deltab-deltaa
      IF(MIN(ABS(res_p),ABS(res_delta_inner)).LT.tol)THEN
        var_m = deltac
        res(1)=res_p
        res(2)=res_delta_inner
        EXIT
      ENDIF
      !check in which interval the zero-crossing takes place
      IF(fa*fc.GT.0.)THEN
        deltaa = deltac
        CALL dT2p_PeTS(pa,deltaa*rhocrit,T)
        fa = pa-p
      ELSEIF(fb*fc.GT.0.) THEN
        deltab = deltac
        CALL dT2p_PeTS(pb,deltab*rhocrit,T)
        fb = pb-p
      ELSE
        WRITE(*,*)'iterate function bisection scheme error in IF(fa*fc.GT.0.)'
        CALL Abort()
      ENDIF
    END DO
  
    delta(1) = var_m
  
    !Calcualte GMASS at liquid and gaseous side
    DO i_inner = 1, 2
      CALL dT2g_PeTS(g(i_inner),delta(i_inner)*rhocrit,T)
    END DO
    !Calculate the resiudal in Gibbs energy
    res_g     = g(2)-g(1)
    res_delta = delta_upper(2)-delta_upper(1)
    !Exit if res is small enough
    IF(ABS(res_delta).LT.tol)THEN
      res=res_delta
      EXIT
    ENDIF
    !set new lower boundary and calculate new DMASSvap(2)
    IF(res_g.LT.0)THEN
      delta_upper(1) = delta(2)
      delta(2) = 0.5*(delta(2)+delta_upper(2))
    !set new upper boundary and calculate new DMASSvap(2)
    ELSE
      delta_upper(2) = delta(2)
      delta(2) = 0.5*(delta(2)+delta_upper(1))
    ENDIF
  ENDDO
ENDIF

! check if iterationnumber is reached
IF (i.GE.EOSIterationMAX .OR. ANY(ABS(res).GT.tol)) THEN
  WRITE(*,*)'IteratePsat failed for T = ',T,'res',res,'i',i
  rho_l = rho_ancillary(1)
  rho_v = rho_ancillary(2)
  CALL AncillaryEquation_Pressure(T,pSat)
  !CALL Abort(!"Saturation conditions not reached in IteratePSat. i,res",i,ABS(res(1))+ABS(res(2)))
ELSE
  ! SI units
  rho_l = delta(1)*rhocrit
  rho_v = delta(2)*rhocrit
  ! calculate pressure
  CALL dT2p_PeTS(pSat,rho_l,T)
ENDIF

END SUBROUTINE IteratePSat


!==================================================================================================================================
!> Iterates the saturation temperature and densities for a given temperature
!==================================================================================================================================
SUBROUTINE IterateTSat(p,TSat,rho_l,rho_v)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! Argument list declaration
! INPUT VARIABLES
REAL,INTENT(IN) :: p
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: TSat
REAL,INTENT(OUT) :: rho_l
REAL,INTENT(OUT) :: rho_v
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION 
INTEGER :: i
REAL    :: pSat
REAL    :: f,res
REAL    :: Ta,Tb
REAL    :: tol
!==================================================================================================================================

IF (p.GE.pcrit) THEN
  rho_l = rhocrit
  rho_v = rho_l
  TSat  = Tcrit
  RETURN
ENDIF

tol = 1.E-8

! Bisection

! Set lower boundary
IF(p.LT.0.9*pcrit)THEN
  Ta = Tmin
ELSE
  Ta = 0.7*Tcrit
ENDIF
! Set critical temperate as upper boundary
Tb = Tcrit

! Start iteration loop
DO i = 1 , EOSIterationMAX
  ! Calculate a new guess
  TSat = 0.5*(Ta+Tb)
  ! Calculate saturation pressure at the guess
  CALL IteratePSat(TSat,pSat,rho_l,rho_v)
  ! Calculate and check residual of the saturation pressure
  f = p - pSat
  res = Ta-Tb
  IF(ABS(res).LT.tol)EXIT
  ! Set new temperature boundaries
  IF(f.GT.0.)THEN
    Ta = TSat
  ELSE
    Tb = TSat
  END IF
END DO

END SUBROUTINE IterateTSat


!==================================================================================================================================!
!==================================================================================================================================!
FUNCTION Eta(tau,delta,ordertau,orderdelta)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: tau
REAL,INTENT(IN)    :: delta
INTEGER,INTENT(IN) :: ordertau
INTEGER,INTENT(IN) :: orderdelta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: Eta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: deriv_pair
REAL    :: d,d_tau,d_tautau
!===================================================================================================================================

! which derivative should be calculated?
deriv_pair = 100*ordertau + orderdelta

! calculate d, if necessary
IF(orderdelta.LT.2) THEN
  d = 1. - c(1)*EXP(-c(2)*tau/Tcrit)
ENDIF

SELECT CASE(deriv_pair)
CASE(0) ! eta
  Eta = pi/6.*delta*rhocrit*d*d*d
CASE(1) ! d(eta)/d(delta)
  Eta = pi/6.*rhocrit*d*d*d
CASE(2) ! d^2(eta)/d(delta)^2
  Eta = 0.
CASE(3) ! d^3(eta)/d(delta)^3
  Eta = 0.
CASE(100) ! d(eta)/d(tau)
  d_tau = c(1)*c(2)/Tcrit*EXP(-c(2)*tau/Tcrit)
  Eta = pi/2.*delta*rhocrit*d*d*d_tau
CASE(200) ! d^2(eta)/d(tau^2)
  d_tau = c(1)*c(2)/Tcrit*EXP(-c(2)*tau/Tcrit)
  d_tautau = -c(1)*c(2)*c(2)/Tcrit/Tcrit*EXP(-c(2)*tau/Tcrit)
  Eta = 0.5*pi*delta*rhocrit*d*(2.*d_tau*d_tau + d*d_tautau)
CASE(101) ! d^2(eta)/d(tau)/d(delta)
  d_tau = c(1)*c(2)/Tcrit*EXP(-c(2)*tau/Tcrit)
  Eta = pi/2.*rhocrit*d*d*d_tau
CASE(102) ! d^3(eta)/d(tau)/d(delta)/d(delta)
  Eta = 0.
CASE DEFAULT
  WRITE(*,*)'Cannot calculate higher order derivative of Eta term. delta=, tau=',orderdelta,ordertau
  CALL Abort()
END SELECT

END FUNCTION Eta



!==================================================================================================================================!
!> Derivatives of the first polynomial term up to arbitrary order
!==================================================================================================================================!
FUNCTION PolyTerm1(tau,delta,order) 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: tau
REAL,INTENT(IN)    :: delta
INTEGER,INTENT(IN) :: order
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: PolyTerm1
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iOrder
REAL    :: coeff
REAL    :: Eta_loc
!===================================================================================================================================

Eta_loc = Eta(tau,delta,0,0)

PolyTerm1 = 0.
DO i = 0, 6
  coeff = 1.
  DO iOrder=1,order
    coeff = coeff*(REAL(i)-(REAL(iOrder)-1.))
  END DO
  PolyTerm1 = PolyTerm1 + coeff*a(i)*Eta_loc**(i-order)
END DO ! i = 0, 6

END FUNCTION PolyTerm1


!==================================================================================================================================!
!> Derivatives of the second polynomial term up to arbitrary order
!==================================================================================================================================!
FUNCTION PolyTerm2(tau,delta,order) 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: tau
REAL,INTENT(IN)    :: delta
INTEGER,INTENT(IN) :: order
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: PolyTerm2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iOrder
REAL    :: coeff
REAL    :: Eta_loc
!===================================================================================================================================

Eta_loc = Eta(tau,delta,0,0)

PolyTerm2 = 0.
DO i = 0, 6
  coeff = 1.
  DO iOrder=1,order
    coeff = coeff*(REAL(i)-(REAL(iOrder)-1.))
  END DO
  PolyTerm2 = PolyTerm2 + coeff*b(i)*Eta_loc**(i-order)
END DO ! i = 0, 6

END FUNCTION PolyTerm2


!==================================================================================================================================!
!> Derivatives of the first term
!==================================================================================================================================!
FUNCTION Term1(tau,delta,order) 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: tau
REAL,INTENT(IN)    :: delta
INTEGER,INTENT(IN) :: order
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: Term1
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: Eta_loc
!===================================================================================================================================

Eta_loc = Eta(tau,delta,0,0)

SELECT CASE(order)
CASE(0)
  Term1 = 3.*Eta_loc/(1.-Eta_loc) + Eta_loc/(1.-Eta_loc)**2
CASE(1)
  Term1 = -2.*(Eta_loc-2.)/(1.-Eta_loc)**3
CASE(2)
  Term1 = (10.-4.*Eta_loc)/(Eta_loc-1.)**4
CASE(3)
  Term1 = 12.*(Eta_loc-3.)/(Eta_loc-1.)**5
CASE DEFAULT
  WRITE(*,*)'Cannot calculate higher order derivative of Term1. order=',order
  CALL Abort()
END SELECT

END FUNCTION Term1


!==================================================================================================================================!
!> Derivatives of the second term
!==================================================================================================================================!
FUNCTION Term2(tau,delta,order) 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: tau
REAL,INTENT(IN)    :: delta
INTEGER,INTENT(IN) :: order
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: Term2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: Eta_loc
!===================================================================================================================================
Eta_loc = Eta(tau,delta,0,0)

SELECT CASE(order)
CASE(0)
  Term2 = (1.+(8*Eta_loc-2.*Eta_loc*Eta_loc)/(1.-Eta_loc)**4)**(-1)
CASE(1)
  Term2 = -4.*(Eta_loc-1.)**3*(Eta_loc**2-5.*Eta_loc-2.)/(Eta_loc**4-4.*Eta_loc**3+4.*Eta_loc**2+4.*Eta_loc+1.)**2
CASE(2)
  Term2 = 4.*(Eta_loc-1.)**2*(3.*Eta_loc**6-30.*Eta_loc**5+77.*Eta_loc**4-80.*Eta_loc**3+39*Eta_loc**2+82*Eta_loc+17.)/ &
          (Eta_loc**4-4.*(Eta_loc**3-Eta_loc**2-Eta_loc)+1.)**3
CASE(3)
  Term2 = -48.*(Eta_loc-1.)*(Eta_loc**10-15.*Eta_loc**9+77.*Eta_loc**8-210.*Eta_loc**7+372.*Eta_loc**6-352.*Eta_loc**5+238*Eta_loc**3-109.*Eta_loc**2-97.*Eta_loc-13.)/  &
          (Eta_loc**4-4.*(Eta_loc**3-Eta_loc**2-Eta_loc)+1.)**4
CASE DEFAULT
  WRITE(*,*)'Cannot calculate higher order derivative of Term2. order=',order
  CALL Abort()
END SELECT

END FUNCTION Term2


!==================================================================================================================================!
!> Return value of the residual free Helmholtz energy
!==================================================================================================================================!
FUNCTION ALPHAR(tau,delta,ordertau,orderdelta) 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: tau
REAL,INTENT(IN)    :: delta
INTEGER,INTENT(IN) :: ordertau
INTEGER,INTENT(IN) :: orderdelta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: ALPHAR
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: deriv_pair
REAL    :: Tc,R
REAL    :: D,D_d
REAL    :: T1,T1_d,T1_dd,T1_ddd
REAL    :: T2,T2_d,T2_dd,T2_ddd
REAL    :: I1,I1_d,I1_dd,I1_ddd
REAL    :: I2,I2_d,I2_dd,I2_ddd
REAL    :: Eta_del,Eta_tau,Eta_tautau,Eta_deltau
REAL    :: aref,a1,a2
REAL    :: TI_d,TI_dd,TI_ddd
!===================================================================================================================================
IF (delta.EQ.0.) THEN
  ALPHAR = 0.
  RETURN
ENDIF

! which derivative should be calculated?
deriv_pair = 100*ordertau + orderdelta

  T1      = Term1(tau,delta,0)
  T1_d    = Term1(tau,delta,1)
  T1_dd   = Term1(tau,delta,2)
  T1_ddd  = Term1(tau,delta,3)
  T2      = Term2(tau,delta,0)
  T2_d    = Term2(tau,delta,1)
  T2_dd   = Term2(tau,delta,2)
  T2_ddd  = Term2(tau,delta,3)
  I1      = PolyTerm1(tau,delta,0)
  I1_d    = PolyTerm1(tau,delta,1)
  I1_dd   = PolyTerm1(tau,delta,2)
  I1_ddd  = PolyTerm1(tau,delta,3)
  I2      = PolyTerm2(tau,delta,0)
  I2_d    = PolyTerm2(tau,delta,1)
  I2_dd   = PolyTerm2(tau,delta,2)
  I2_ddd  = PolyTerm2(tau,delta,3)
  D       = pi*delta*rhocrit
  D_d     = pi*rhocrit
  Eta_del = Eta(tau,delta,0,1)
  Eta_tau = Eta(tau,delta,1,0)
  Eta_tautau = Eta(tau,delta,2,0)
  Eta_deltau = Eta(tau,delta,1,1)

  Tc = Tcrit
  R  = R_mass

SELECT CASE(deriv_pair)
CASE(0) ! alphar
  aref = T1/R
  a1   = -2./(R*Tc)*tau*D*I1
  a2   = -1./(R*Tc*Tc)*tau*tau*D*T2*I2
CASE(1) ! d(alphar)/d(delta)
  aref = T1_d*Eta_del/R
  a1   = -2.*tau/(R*Tc)*(D_d*I1 + D*I1_d*Eta_del)
  TI_d = T2_d*I2+T2*I2_d
  a2   = -tau*tau/(R*Tc*Tc)*(D_d*T2*I2 + D*Eta_del*TI_d)
CASE(2) ! d^2(alphar)/d(delta)^2
  TI_d  = T2_d*I2+T2*I2_d
  TI_dd = T2_dd*I2 + 2.*T2_d*I2_d + T2*I2_dd
  aref = T1_dd*Eta_del*Eta_del/R
  a1   = -2.*tau/(R*Tc)*Eta_del*(2.*D_d*I1_d + D*I1_dd*Eta_del)
  a2   = -tau*tau/(R*Tc*Tc)*Eta_del*(2.*D_d*TI_d + D*Eta_del*TI_dd)
CASE(3) ! d^3(alphar)/d(delta)^3
  TI_dd  = T2_dd*I2 + 2.*T2_d*I2_d + T2*I2_dd
  TI_ddd = T2_ddd*I2 + 3.*(T2_dd*I2_d + T2_d*I2_dd) + T2*I2_ddd
  aref = T1_ddd*Eta_del*Eta_del*Eta_del/R
  a1 = -2.*tau/(R*Tc)*Eta_del*Eta_del*(3.*D_d*I1_dd + D*I1_ddd*Eta_del)
  a2 = -tau*tau/(R*Tc*Tc)*Eta_del*Eta_del*(3.*D_d*TI_dd + D*Eta_del*TI_ddd)
CASE(100) ! d(alphar)/d(tau)
  TI_d  = T2_d*I2+T2*I2_d
  aref = T1_d*Eta_tau/R
  a1 = -2.*D/(R*Tc)*(I1 + tau*I1_d*Eta_tau)
  a2 = -D/(R*Tc*Tc)*tau*(2.*T2*I2 + tau*Eta_tau*TI_d)
CASE(200) ! d^2(alphar)/d(tau^2)
  TI_d  = T2_d*I2+T2*I2_d
  TI_dd = T2_dd*I2 + 2.*T2_d*I2_d + T2*I2_dd
  aref = 1./R*(T1_dd*Eta_tau*Eta_tau + T1_d*Eta_tautau)
  a1 = -2.*D/(R*Tc)*(2.*I1_d*Eta_tau + tau*(I1_dd*Eta_tau*Eta_tau + I1_d*Eta_tautau))
  a2 = -D/(R*Tc*Tc)*(2*T2*I2 + TI_d*tau*(4.*Eta_tau + tau*Eta_tautau) + tau*tau*Eta_tau*Eta_tau*TI_dd)
CASE(101) ! d^2(alphar)/d(tau)/d(delta)
  TI_d  = T2_d*I2+T2*I2_d
  TI_dd = T2_dd*I2 + 2.*T2_d*I2_d + T2*I2_dd
  aref = (T1_dd*Eta_del*Eta_tau + T1_d*Eta_deltau)/R
  a1 = -2./(R*Tc)*((D_d*I1 + D*I1_d*Eta_del) + tau*(D_d*I1_d*Eta_tau + D*(I1_dd*Eta_del*Eta_tau + I1_d*Eta_deltau)))
  a2 = -tau/(R*Tc*Tc)*(2.*D_d*T2*I2 + TI_d*(2.*D*Eta_del + tau*(D_d*Eta_tau + D*Eta_deltau)) + D*tau*Eta_del*Eta_tau*TI_dd)
CASE(102) ! d^3(alphar)/d(tau)/d(delta)/d(delta)
  TI_d  = T2_d*I2+T2*I2_d
  TI_dd = T2_dd*I2 + 2.*T2_d*I2_d + T2*I2_dd
  TI_ddd = T2_ddd*I2 + 3.*(T2_dd*I2_d + T2_d*I2_dd) + T2*I2_dd
  aref = (T1_ddd*Eta_del*Eta_del*Eta_tau + T1_dd*2.*Eta_del*Eta_deltau)/R
  a1 = -2./(R*Tc)*((2.*D_d*I1_d + D*I1_dd*Eta_del)*(Eta_del + tau*Eta_deltau) + &
                   tau*Eta_del*(2.*D_d*I1_dd*Eta_tau + D*(I1_ddd*Eta_del*Eta_tau + I1_dd*Eta_deltau)))
  a2 = -1./(R*Tc*Tc)*(2.*D_d*TI_d + D*Eta_del*TI_dd*(2.*tau*Eta_del + tau*tau*Eta_deltau) + &
                      tau*tau*Eta_del*(2.*D_d*Eta_tau*TI_dd + D*Eta_deltau*TI_dd + D*Eta_del*TI_ddd))
CASE DEFAULT
  WRITE(*,*)'Cannot calculate higher order derivative of alphar term. delta=, tau=',orderdelta,ordertau
  CALL Abort()
END SELECT

ALPHAR = aref + a1 + a2
END FUNCTION ALPHAR


!==================================================================================================================================!
!> Return value of the ideal gas contribution of the free Helmholtz energy
!>
!> The function for the ideal gas contribution is
!> alpha0 = ln(delta) + c(tau) + f(h0,s0)
!==================================================================================================================================!
FUNCTION ALPHA0(tau,delta,ordertau,orderdelta) 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: tau
REAL,INTENT(IN)    :: delta
INTEGER,INTENT(IN) :: ordertau
INTEGER,INTENT(IN) :: orderdelta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: ALPHA0
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: deriv_pair
REAL    :: delta0,tau0
!===================================================================================================================================

! which derivative should be calculated?
deriv_pair = 100*ordertau + orderdelta
tau0   = Tcrit/T0
delta0 = rho0/rhocrit

SELECT CASE(deriv_pair)
CASE(0) ! alpha0
  IF (delta.EQ.0.) THEN
    ALPHA0 = -1.E12
  ELSE
    ALPHA0 = LOG(delta) + 1.5*LOG(tau) + ig1*tau + ig2
  ENDIF
CASE(1) ! d(alpha0)/d(delta)
  IF (delta.EQ.0.) THEN
    ALPHA0 = 1.E12
  ELSE
    ALPHA0 = 1./delta
  ENDIF
CASE(2) ! d^2(alpha0)/d(delta)^2
  IF (delta.EQ.0.) THEN
    ALPHA0 = -1.E12
  ELSE
    ALPHA0 = -1./(delta*delta)
  ENDIF
CASE(100) ! d(alpha0)/d(tau)
  ALPHA0 = 1.5/tau + ig1
CASE(200) ! d^2(alpha0)/d(tau^2)
  ALPHA0 = -1.5/(tau*tau)
CASE(101) ! d^2(alpha0)/d(tau)/d(delta)
  ALPHA0 = 0.
CASE DEFAULT
  WRITE(*,*)'No such derivative of alpha_i can be calculated',deriv_pair
  CALL Abort()
END SELECT

END FUNCTION ALPHA0


!==================================================================================================================================!
!> (rho,T) -> (u)
!==================================================================================================================================!
SUBROUTINE dT2u_PeTS(u,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: u
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphai10
REAL :: alphar10
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
!IF(delta.LT.0. .OR. tau.LT.0. .OR. T.LT.Tmin(1)) THEN
  !RETURN
!ENDIF
IF(T.EQ.0.) THEN
  u = Tcrit*R_mass*ig1
  RETURN
ENDIF
! internal energy
alphai10 = ALPHA0(tau,delta,1,0)
alphar10 = ALPHAR(tau,delta,1,0)
u        = R_mass*T*tau*(alphai10+alphar10)
END SUBROUTINE dT2u_PeTS


!==================================================================================================================================!
!> (rho,T) -> (p)
!==================================================================================================================================!
SUBROUTINE dT2p_PeTS(p,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: p
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphar01
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
IF(T.EQ.0.) THEN
  p = 0.0
  RETURN
ENDIF
! pressure
alphar01 = ALPHAR(tau,delta,0,1)
p        = rho*T*R_mass*(1.0 + delta*alphar01)
END SUBROUTINE dT2p_PeTS


!==================================================================================================================================!
!> (rho,T) -> (h)
!==================================================================================================================================!
SUBROUTINE dT2h_PeTS(h,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: h
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphai10
REAL :: alphar10,alphar01
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
!IF(delta.LT.0. .OR. tau.LT.0. .OR. T.LT.Tmin(1)) THEN
  !RETURN
!ENDIF
IF(T.EQ.0.) THEN
  h = Tcrit*R_mass*ig1
  RETURN
ENDIF
! enthalpy
alphai10 = ALPHA0(tau,delta,1,0)
alphar10 = ALPHAR(tau,delta,1,0)
alphar01 = ALPHAR(tau,delta,0,1)
h        = T*R_mass*((1.0 + delta*alphar01) + tau*(alphai10+alphar10))
END SUBROUTINE dT2h_PeTS


!==================================================================================================================================!
!> (rho,T) -> (a)
!==================================================================================================================================!
SUBROUTINE dT2a_PeTS(a,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: a
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphai20
REAL :: alphar01,alphar02,alphar11,alphar20
REAL :: dpdrho_T,dpdT_rho,dsdrho_T,dsdT_rho
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
!IF(delta.LT.0. .OR. tau.LT.0. .OR. T.LT.Tmin(1)) THEN
  !a = 0.
  !RETURN
!ENDIF
alphai20 = ALPHA0(tau,delta,2,0)
alphar20 = ALPHAR(tau,delta,2,0)
dsdT_rho = R_mass/T*(-tau*tau*(alphai20 + alphar20))
IF (delta.EQ.0.) THEN
  a        = R_mass*(T+R_mass/dsdT_rho)
ELSE
  !speed of sound
  alphar01 = ALPHAR(tau,delta,0,1)
  alphar02 = ALPHAR(tau,delta,0,2)
  dpdrho_T = T*R_mass*(1. + 2.*delta*alphar01 + delta*delta*alphar02)
  alphar11 = ALPHAR(tau,delta,1,1)
  dpdT_rho = rho*R_mass*(1. + delta*alphar01 - tau*delta*alphar11)
  dsdrho_T = R_mass/rho*(-1.*(1. + delta*alphar01 - tau*delta*alphar11))
  a        = dpdrho_T - dpdT_rho*dsdrho_T/dsdT_rho
ENDIF
IF (a.GT.0.) THEN
  a = SQRT(a)
ELSE
  a = 0.
ENDIF
END SUBROUTINE dT2a_PeTS


!==================================================================================================================================!
!> (rho,T) -> (cp)
!==================================================================================================================================!
SUBROUTINE dT2cp_PeTS(cp,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: cp
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphai20
REAL :: alphar01,alphar02,alphar11,alphar20
REAL :: dhdT_rho,dhdrho_T
REAL :: dpdrho_T,dpdT_rho
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
!IF(delta.LT.0. .OR. tau.LT.0. .OR. T.LT.Tmin(1)) THEN
  !RETURN
!ENDIF
IF(T.EQ.0.) THEN
  cp = 1.5*R_mass
  RETURN
ENDIF
! cp
alphai20 = ALPHA0(tau,delta,2,0)
alphar20 = ALPHAR(tau,delta,2,0)
alphar01 = ALPHAR(tau,delta,0,1)
alphar11 = ALPHAR(tau,delta,1,1)
dhdT_rho = R_mass*(-tau*tau*(alphai20 + alphar20) + (1. + delta*alphar01 - tau*delta*alphar11))
alphar02 = ALPHAR(tau,delta,0,2)
IF(delta.EQ.0.) THEN
  dhdrho_T = 0.
ELSE
  dhdrho_T = T*R_mass/rho*(tau*delta*alphar11 + delta*alphar01 + delta*delta*alphar02)
ENDIF
dpdT_rho = rho*R_mass*(1. + delta*alphar01 - tau*delta*alphar11)
dpdrho_T = T*R_mass*(1. + 2.*delta*alphar01 + delta*delta*alphar02)
cp       = dhdT_rho - dhdrho_T*dpdT_rho/dpdrho_T
END SUBROUTINE dT2cp_PeTS


!==================================================================================================================================!
!> (rho,T) -> (cv)
!==================================================================================================================================!
SUBROUTINE dT2cv_PeTS(cv,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: cv
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphai20
REAL :: alphar20
REAL :: dudT_rho
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
!IF(delta.LT.0. .OR. tau.LT.0. .OR. T.LT.Tmin(1)) THEN
  !RETURN
!ENDIF
IF(T.EQ.0.) THEN
  cv = 1.5*R_mass
  RETURN
ENDIF
! cv
alphai20 = ALPHA0(tau,delta,2,0)
alphar20 = ALPHAR(tau,delta,2,0)
dudT_rho = R_mass*(-tau*tau*(alphai20 + alphar20))
cv       = dudT_rho
END SUBROUTINE dT2cv_PeTS


!==================================================================================================================================!
!> (rho,T) -> (s)
!==================================================================================================================================!
SUBROUTINE dT2s_PeTS(s,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: s
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphai00,alphai10
REAL :: alphar00,alphar10
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
!IF(delta.LT.0. .OR. tau.LT.0.) THEN
  !!s = 0.
  !RETURN
!ENDIF
IF(T.EQ.0.) THEN
  s = 0.
  RETURN
ENDIF
! entropy
alphai00 = ALPHA0(tau,delta,0,0)
alphai10 = ALPHA0(tau,delta,1,0)
alphar00 = ALPHAR(tau,delta,0,0)
alphar10 = ALPHAR(tau,delta,1,0)
s        = R_mass*(tau*(alphai10+alphar10) - (alphai00+alphar00))
END SUBROUTINE dT2s_PeTS


!==================================================================================================================================!
!> (rho,T) -> (g)
!==================================================================================================================================!
SUBROUTINE dT2g_PeTS(g,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: g
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphai00
REAL :: alphar00,alphar01
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
!IF(delta.LT.0. .OR. tau.LT.0. .OR. T.LT.Tmin(1)) THEN
  !RETURN
!ENDIF
IF(T.EQ.0.) THEN
  g = Tcrit*R_mass*ig1
  RETURN
ENDIF
! gibbs energy
alphai00 = ALPHA0(tau,delta,0,0)
alphar00 = ALPHAR(tau,delta,0,0)
alphar01 = ALPHAR(tau,delta,0,1)
g        = T*R_mass*((1.0 + delta*alphar01) + (alphai00+alphar00))
END SUBROUTINE dT2g_PeTS


!==================================================================================================================================!
!> (rho,T) -> (phase identifier)
!> formula found in
!> Venkatarathnam and Oellrich, Fluid Phase Equilibria 2011
!> https://doi.org/10.1016/j.fluid.2010.12.001
!==================================================================================================================================!
SUBROUTINE dT2pip_PeTS(pip,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT) :: pip
REAL,INTENT(IN)  :: rho
REAL,INTENT(IN)  :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphar01,alphar02
REAL :: alphar11,alphar12
REAL :: d2pdrhodT,dpdT_rho,d2pdrho2_T,dpdrho_T
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
!IF(delta.LT.0. .OR. tau.LT.0.) THEN
  !pip = -999.
  !RETURN
!ENDIF
IF(delta.EQ.0. .OR. T.EQ.0.) THEN
  pip = 0.
  RETURN
ENDIF
IF (T.LT.Tmin) THEN
  tau = Tcrit/Tmin
ENDIF
! phase identification parameter
alphar01 = ALPHAR(tau,delta,0,1)
alphar02 = ALPHAR(tau,delta,0,2)
alphar11 = ALPHAR(tau,delta,1,1)
alphar12 = ALPHAR(tau,delta,1,2)
CALL dpdT_rho_PeTS(dpdT_rho,rho,T)
CALL dpdrho_T_PeTS(dpdrho_T,rho,T)
CALL d2pdrho2_T_PeTS(d2pdrho2_T,rho,T)
d2pdrhodT = R_mass*(1. + 2.*delta*alphar01 + delta*delta*alphar02 - 2.*delta*tau*alphar11 - tau*delta*delta*alphar12)
pip = 2. - rho*(d2pdrhodT/dpdT_rho - d2pdrho2_T/dpdrho_T)
END SUBROUTINE dT2pip_PeTS


!==================================================================================================================================!
!> derivative d(p)/d(T) at constant rho
!==================================================================================================================================!
SUBROUTINE dpdT_rho_PeTS(dpdT_rho,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: dpdT_rho
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphar11,alphar01
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
!IF(delta.LT.0. .OR. tau.LT.0. .OR. T.LT.Tmin(1)) THEN
  !RETURN
!ENDIF
! pressure derivative
alphar01 = ALPHAR(tau,delta,0,1)
alphar11 = ALPHAR(tau,delta,1,1)
dpdT_rho = rho*R_mass*(1. + delta*alphar01 - tau*delta*alphar11)
END SUBROUTINE dpdT_rho_PeTS


!==================================================================================================================================!
!> derivative d(u)/d(T) at constant rho
!==================================================================================================================================!
SUBROUTINE dudT_rho_PeTS(dudT_rho,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: dudT_rho
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL dT2cv_PeTS(dudT_rho,rho,T)
END SUBROUTINE dudT_rho_PeTS

!==================================================================================================================================!
!> derivative d(u)/d(rho) at constant T
!==================================================================================================================================!
SUBROUTINE dudrho_T_PeTS(dudT_rho,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: dudT_rho
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphai20
REAL :: alphar20
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
!IF(delta.LT.0. .OR. tau.LT.0. .OR. T.LT.Tmin(1)) THEN
  !RETURN
!ENDIF
alphai20 = ALPHA0(tau,delta,2,0)
alphar20 = ALPHAR(tau,delta,2,0)
dudT_rho = R_mass*(-tau*tau*(alphai20 + alphar20))
END SUBROUTINE dudrho_T_PeTS

!==================================================================================================================================!
!> derivative d(h)/d(T) at constant rho
!==================================================================================================================================!
SUBROUTINE dhdT_rho_PeTS(dhdT_rho,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: dhdT_rho
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphai20
REAL :: alphar01,alphar11,alphar20
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
!IF(delta.LT.0. .OR. tau.LT.0. .OR. T.LT.Tmin(1)) THEN
  !RETURN
!ENDIF
alphai20 = ALPHA0(tau,delta,2,0)
alphar20 = ALPHAR(tau,delta,2,0)
alphar01 = ALPHAR(tau,delta,0,1)
alphar11 = ALPHAR(tau,delta,1,1)
dhdT_rho = R_mass*(-tau*tau*(alphai20 + alphar20) + (1. + delta*alphar01 - tau*delta*alphar11))
END SUBROUTINE dhdT_rho_PeTS


!==================================================================================================================================!
!> derivative d(s)/d(T) at constant rho
!==================================================================================================================================!
SUBROUTINE dsdT_rho_PeTS(dsdT_rho,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: dsdT_rho
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphai20
REAL :: alphar20
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
!IF(delta.LT.0. .OR. tau.LT.0. .OR. T.LT.Tmin(1)) THEN
  !RETURN
!ENDIF
alphai20 = ALPHA0(tau,delta,2,0)
alphar20 = ALPHAR(tau,delta,2,0)
dsdT_rho = R_mass/T*(-tau*tau*(alphai20 + alphar20))
END SUBROUTINE dsdT_rho_PeTS


!==================================================================================================================================!
!> derivative d(p)/d(rho) at constant T
!==================================================================================================================================!
SUBROUTINE dpdrho_T_PeTS(dpdrho_T,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: dpdrho_T
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphar01,alphar02
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
!IF(delta.LT.0. .OR. tau.LT.0. .OR. T.LT.Tmin(1)) THEN
  !RETURN
!ENDIF
! pressure derivative
alphar01 = ALPHAR(tau,delta,0,1)
alphar02 = ALPHAR(tau,delta,0,2)
dpdrho_T = T*R_mass*(1. + 2.*delta*alphar01 + delta*delta*alphar02)
END SUBROUTINE dpdrho_T_PeTS


!==================================================================================================================================!
!> second derivative d^2(p)/d(rho^2) at constant T
!==================================================================================================================================!
SUBROUTINE d2pdrho2_T_PeTS(d2pdrho2_T,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: d2pdrho2_T
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: tau,delta
REAL :: alphar01,alphar02,alphar03
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
IF(delta.EQ.0.) THEN
  d2pdrho2_T=0.
  RETURN
ENDIF
! pressure derivative
alphar01 = ALPHAR(tau,delta,0,1)
alphar02 = ALPHAR(tau,delta,0,2)
alphar03 = ALPHAR(tau,delta,0,3)
d2pdrho2_T = T*R_mass/rho*(2.*delta*alphar01 + 4.*delta*delta*alphar02 + delta*delta*delta*alphar03)
END SUBROUTINE d2pdrho2_T_PeTS


! Derivatives, that cannot be calculated directly, but with chain rule etc. 
!==================================================================================================================================!
!> derivative d(p)/d(rho) at constant u
!==================================================================================================================================!
SUBROUTINE dpdrho_u_PeTS(dpdrho_u,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: dpdrho_u
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                         :: tau,delta
REAL                                         :: dpdrho_T,dpdT_rho,dudrho_T,dudT_rho
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
! Calc needed derivatives
CALL dpdrho_T_PeTS(dpdrho_T,rho,T)
CALL dpdT_rho_PeTS(dpdT_rho,rho,T)
CALL dudrho_T_PeTS(dudrho_T,rho,T)
CALL dudT_rho_PeTS(dudT_rho,rho,T)
! chain rule
dpdrho_u = dpdrho_T - dpdT_rho * dudrho_T * (1./dudT_rho)
END SUBROUTINE dpdrho_u_PeTS
!==================================================================================================================================!
!> derivative d(p)/d(u) at constant rho
!==================================================================================================================================!
SUBROUTINE dpdu_rho_PeTS(dpdu_rho,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: dpdu_rho
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                         :: tau,delta
REAL                                         :: dpdT_rho,dudT_rho
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
! Calc needed derivatives
CALL dpdT_rho_PeTS(dpdT_rho,rho,T)
CALL dudT_rho_PeTS(dudT_rho,rho,T)
! chain rule
dpdu_rho = dpdT_rho * (1./dudT_rho) 
END SUBROUTINE dpdu_rho_PeTS
!==================================================================================================================================!
!> derivative d(T)/d(rho) at constant u
!==================================================================================================================================!
SUBROUTINE dTdrho_u_PeTS(dTdrho_u,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: dTdrho_u
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                         :: tau,delta
REAL                                         :: dudrho_T,dudT_rho
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
! Calc needed derivatives
CALL dudrho_T_PeTS(dudrho_T,rho,T)
CALL dudT_rho_PeTS(dudT_rho,rho,T)
! chain rule
dTdrho_u = -dudrho_T * (1./dudT_rho) 
END SUBROUTINE dTdrho_u_PeTS
!==================================================================================================================================!
!> derivative d(T)/d(u) at constant rho
!==================================================================================================================================!
SUBROUTINE dTdu_rho_PeTS(dTdu_rho,rho,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: dTdu_rho
REAL,INTENT(IN)                              :: rho
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                         :: tau,delta
REAL                                         :: dudT_rho
!===================================================================================================================================
! calculate reduced variables
tau   = Tcrit/T
delta = rho/rhocrit
! Calc needed derivatives
CALL dudT_rho_PeTS(dudT_rho,rho,T)
! chain rule
dTdu_rho = (1./dudT_rho) 
END SUBROUTINE dTdu_rho_PeTS


!==================================================================================================================================!
!> calculate the spinodal densities

!> find d(p)/d(rho)=0 at constant T
!> Note, for a general EOS (or a multiparameter EOS), the only valid spinodals are the first ones coming from the vapour and 
!> the liquid sides!
!==================================================================================================================================!
SUBROUTINE CalcSpinodal_PeTS(rho_l,rho_v,T)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(OUT)                             :: rho_l
REAL,INTENT(OUT)                             :: rho_v
REAL,INTENT(IN)                              :: T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,i_inner
REAL    :: rho,drho,rhoc,delta,f,res,dvar_iterdrho,var_iter
REAL    :: tol
!===================================================================================================================================
! find d(p)/d(rho)=0 at constant T

! set local tolerance
tol = EOSIterationTOL

! vapour spinodal
rho = 0.001 ! choose a very small density
DO i = 1, EOSIterationMAX
  IF (i.EQ.10) tol=10.*tol
  IF (i.EQ.30) tol=10.*tol
  IF (i.EQ.50) tol=10.*tol
  CALL dpdrho_T_PeTS(var_iter,rho,T) ! target function
  f = var_iter
  res  =  ABS(f)
  IF(res.LT.tol) EXIT
  CALL d2pdrho2_T_PeTS(dvar_iterdrho,rho,T) ! derivative of target function
  drho = -f/dvar_iterdrho
  delta= 0.9
  ! Globalization strategy: linesearch   
  DO i_inner=1,10
    rhoc = rho+delta*drho
    CALL dpdrho_T_PeTS(var_iter,rhoc,T)
    IF((ABS(var_iter).LT.res)) THEN
      EXIT
    ELSE
      delta = 0.5*delta
    ENDIF
  ENDDO
  rho=rho+delta*drho
ENDDO
rho_v = rho

! liquid spinodal
! set local tolerance
tol = EOSIterationTOL
rho = rhomax ! choose a very large density
DO i = 1, EOSIterationMAX
  IF (i.EQ.10) tol=10.*tol
  IF (i.EQ.30) tol=10.*tol
  IF (i.EQ.50) tol=10.*tol
  CALL dpdrho_T_PeTS(var_iter,rho,T) ! target function
  f = var_iter
  res  =  ABS(f)
  IF(res.LT.tol) EXIT
  CALL d2pdrho2_T_PeTS(dvar_iterdrho,rho,T) ! derivative of target function
  drho = -f/dvar_iterdrho
  delta= 0.9
  ! Globalization strategy: linesearch   
  DO i_inner=1,10
    rhoc = rho+delta*drho
    CALL dpdrho_T_PeTS(var_iter,rhoc,T) ! target function
    IF((ABS(var_iter).LT.res)) THEN
      EXIT
    ELSE
      delta = 0.5*delta
    ENDIF
  ENDDO
  rho=rho+delta*drho
ENDDO
rho_l = rho

END SUBROUTINE CalcSpinodal_PeTS


!==================================================================================================================================!
!> Initial guess for temperature
!==================================================================================================================================!
SUBROUTINE InitialGuess_Default_PeTS(T_min,T_max,T,rho,p)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)  :: rho
REAL,INTENT(IN)  :: p
REAL,INTENT(OUT) :: T_min,T_max,T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
T_min = Tmin
T_max = Tmax
T = 0.5*(T_min+T_max)

END SUBROUTINE InitialGuess_Default_PeTS


!==================================================================================================================================!
!> Initial guess for temperature using the Van-Der-Waal law
!> when \rho and p are known
!==================================================================================================================================!
SUBROUTINE InitialGuess_VdW_PeTS(T_min,T_max,T,rho,p)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)  :: rho
REAL,INTENT(IN)  :: p
REAL,INTENT(OUT) :: T_min,T_max,T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: a_VdW,b_VdW
!===================================================================================================================================
T_min = Tmin
T_max = Tmax

! VdW constants
a_VdW = 3.*pcrit/(rhocrit*rhocrit)
b_VdW = 1./(3.*rhocrit)

IF((p.GT.0.) .AND. (rho.LT.rhocrit)) THEN
  T = (p+a_VdW*rhocrit*rhocrit)*(1./rhocrit-b_VdW)/R_mass
ELSE
  T = Tcrit
ENDIF

END SUBROUTINE InitialGuess_VdW_PeTS


!==================================================================================================================================!
!> Initial guess for temperature if internal energy is given
!==================================================================================================================================!
SUBROUTINE InitialGuess_InternalEnergy_PeTS(T_min,T_max,T,rho,u)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)  :: rho
REAL,INTENT(IN)  :: u
REAL,INTENT(OUT) :: T_min,T_max,T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: Tcrit,Tmean
REAL :: ucrit
!===================================================================================================================================
! get some temperatures
T_min = Tmin
T_max = Tmax
Tmean = 0.5*(T_min+T_max)
Tcrit = Tcrit

! check energy for guess
CALL dT2u_pets(ucrit,rho,Tcrit)

! assume a monotone function u(T)
IF (u.GT.ucrit) THEN
  T_min = Tcrit
  T_max = T_max
  T = 0.5*(T_min+T_max)
ELSE
  T_min = T_min
  T_max = Tcrit
  T = 0.5*(T_min+T_max)
ENDIF

END SUBROUTINE InitialGuess_InternalEnergy_PeTS


!==================================================================================================================================!
!> Initial guess for temperature if entropy is given
!==================================================================================================================================!
SUBROUTINE InitialGuess_Entropy_PeTS(T_min,T_max,T,rho,s)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)  :: rho
REAL,INTENT(IN)  :: s
REAL,INTENT(OUT) :: T_min,T_max,T
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: Tcrit,Tmean
REAL :: scrit
!===================================================================================================================================
! get some temperatures
T_min = Tmin
T_max = Tmax
Tmean = 0.5*(T_min+T_max)
Tcrit = Tcrit

! check energy for guess
CALL dT2s_pets(scrit,rho,Tcrit)

! assume a monotone function s(T)
IF (s.GT.scrit) THEN
  T_min = Tcrit
  T_max = T_max
  T = 0.5*(T_min+T_max)
ELSE
  T_min = T_min
  T_max = Tcrit
  T = 0.5*(T_min+T_max)
ENDIF

END SUBROUTINE InitialGuess_Entropy_PeTS


!==================================================================================================================================!
!> Ancillary equation to calculate the saturated density of the vapour and liquid phases
!> 
!> https://doi.org/10.1007/s10765-014-1764-4
!> eqs. (11) and (12)
!> ATTN: The equations were interchanged in the paper.
!> The coefficients were consequently fitted to the wrong equations.
!>
!> https://doi.org/10.1063/1.4945000
!> eqs. (45) - (46)
!==================================================================================================================================!
SUBROUTINE AncillaryEquation_DensitySat(T,rho_sat)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: T
REAL,INTENT(OUT)   :: rho_sat(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
REAL    :: theta
REAL    :: rho_vap,rho_liq
!===================================================================================================================================
theta = (1.-T/Tcrit)

rho_liq = 0.
DO i=1,4
  rho_liq = rho_liq + N_liq(i)*theta**t_liq(i)
END DO
rho_liq = EXP(rho_liq)*rhocrit

rho_vap = 0.
DO i=1,4
  rho_vap = rho_vap + N_vap(i)*theta**t_vap(i)
END DO
rho_vap = (1.+rho_vap)*rhocrit

rho_sat(1) = rho_liq
rho_sat(2) = rho_vap

END SUBROUTINE AncillaryEquation_DensitySat


!==================================================================================================================================!
!> Ancillary equation to calculate the saturated pressure
!> Modified Wagner Equation
!> 
!> https://doi.org/10.1007/s10765-014-1764-4
!> eqs. (10)
!>
!> https://doi.org/10.1063/1.4945000
!> eqs. (44)
!==================================================================================================================================!
SUBROUTINE AncillaryEquation_Pressure(T,p_sat)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: T
REAL,INTENT(OUT)   :: p_sat
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
REAL    :: theta
!===================================================================================================================================
theta = (1.-T/Tcrit)

p_sat = 0.
DO i=1,5
  p_sat = p_sat + N_pres(i)*theta**t_a(i)
END DO
p_sat = p_sat*Tcrit/T
p_sat = EXP(p_sat)*pcrit

END SUBROUTINE AncillaryEquation_Pressure

END MODULE PETS
