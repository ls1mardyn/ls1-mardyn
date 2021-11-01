#include "pets.h"

!==================================================================================================================================
!> Description
!==================================================================================================================================
SUBROUTINE PETSEOS(n_z,integer_x,var_x,integer_y,var_y,integer_z,var_z)
! MODULES                                                                                                                          
  use pets, ONLY: Init,Evaluate
!----------------------------------------------------------------------------------------------------------------------------------
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,INTENT(IN)  :: n_z
INTEGER,INTENT(IN)  :: integer_x,integer_y
INTEGER,INTENT(IN)  :: integer_z(n_z)
DOUBLE PRECISION,INTENT(IN)     :: var_x,var_y
DOUBLE PRECISION,INTENT(OUT)    :: var_z(n_Z)
!===================================================================================================================================
CALL Init()

CALL Evaluate(n_z,integer_x,var_x,integer_y,var_y,integer_z,var_z)
END SUBROUTINE PETSEOS
