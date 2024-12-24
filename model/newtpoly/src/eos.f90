SUBROUTINE EOS_PRESSURE (j_in, k_in, l_in, rho_in, eps_in, p_out)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

INTEGER, INTENT (IN) :: j_in, k_in, l_in
REAL*8, INTENT (IN) :: rho_in, eps_in
REAL*8, INTENT (OUT) :: p_out

p_out = rho_in*eps_in*(ggas-1.0d0)

END SUBROUTINE

SUBROUTINE EOS_SOUNDSPEED (p_in, rho_in, eps_in, cs_out)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

REAL*8, INTENT (IN) :: p_in, rho_in, eps_in
REAL*8, INTENT (OUT) :: cs_out

cs_out = DSQRT(ggas*p_in/rho_in)

END SUBROUTINE

SUBROUTINE EOS_EPSILON (p_in, rho_in, eps_out)
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

REAL*8, INTENT (IN) :: p_in, rho_in
REAL*8, INTENT (OUT) :: eps_out

eps_out = p_in/rho_in/(ggas - 1.0D0)

END SUBROUTINE