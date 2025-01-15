SUBROUTINE EOS_PRESSURE (j_in, k_in, l_in, rho_in, eps_in, p_out)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

INTEGER, INTENT (IN) :: j_in, k_in, l_in
REAL*8, INTENT (IN) :: rho_in, eps_in
REAL*8, INTENT (OUT) :: p_out

p_out = kpoly*rho_in**ggas
CALL EOS_EPSILON(p_out, rho_in, eps(j_in, k_in, l_in))

END SUBROUTINE

SUBROUTINE EOS_SOUNDSPEED (p_in, rho_in, eps_in, cs_out)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

REAL*8, INTENT (IN) :: p_in, rho_in, eps_in
REAL*8, INTENT (OUT) :: cs_out

cs_out = DSQRT(ggas*kpoly*rho_in)

END SUBROUTINE

SUBROUTINE EOS_EPSILON (p_in, rho_in, eps_out)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

REAL*8, INTENT (IN) :: p_in, rho_in
REAL*8, INTENT (OUT) :: eps_out

eps_out = kpoly*rho_in

END SUBROUTINE