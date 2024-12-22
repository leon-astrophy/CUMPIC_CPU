!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine contains the essential for the calculation of      !
! pressure, sound speed and other thermodynamic quantities.          !
! The EOSEPSILON and EOSSOUNDSPEED are required in reconstruction    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOS_PRESSURE (j_in, k_in, l_in, rho_in, eps_in, p_out)
!$ACC ROUTINE SEQ
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

! Input !
INTEGER, INTENT (IN) :: j_in, k_in, l_in

! Input density !
REAL*8, INTENT (IN) :: rho_in, eps_in

! Output value !
REAL*8, INTENT (OUT) :: p_out

REAL*8 :: xe

!-----------------------------------------------------!

! find pressure !
xe = (rho_in/bmax)**(1.0D0/3.0D0)
IF (xe > 1.0D-2) THEN
	p_out = amax*large_pressure(xe)
ELSE
	p_out = amax*small_pressure(xe)
END IF
!------------------------------------------------------!
contains

  REAL*8 function large_pressure(x)
  !$acc routine seq
	implicit none
	REAL*8 :: x
	large_pressure = x*DSQRT(x**2 + 1.0D0)*(2.0D0*x**2 - 3.0D0) + 3.0D0*log(x + DSQRT(x**2 + 1.0D0))
	end function

	REAL*8 function small_pressure(x)
  !$acc routine seq
	implicit none
	REAL*8 :: x
	small_pressure = 1.6D0*x**5 - (4.0D0/7.0D0)*x**7 + (1.0D0/3.0D0)*x**9 - (5.0D0/2.2D1)*x**11 & 
			+ (3.5D1/2.08D2)*x**13 - (2.1D1/1.6D2)*x**15 + (2.31D2/2.176D3)*x**17
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EOS_EPSILON (rho_in, p_in, eps_out)
!$acc routine seq
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

REAL*8 :: xe
REAL*8, INTENT (IN) :: rho_in, p_in
REAL*8, INTENT (OUT) :: eps_out

xe = (rho_in/bmax)**(1.0D0/3.0D0)
IF(xe > 1.0D-2) THEN
	eps_out = amax*large_energy(xe)/rho_in
ELSE
	eps_out = amax*small_energy(xe)/rho_in
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	REAL*8 function large_energy(x)
	!$acc routine seq
	implicit none
	REAL*8 :: x
	large_energy = 3.0D0*x*DSQRT(x**2 + 1.0D0)*(1.0D0 + 2.0D0*x**2) - 3.0D0*log(x + DSQRT(x**2 + 1.0D0))
	end function

	REAL*8 function small_energy(x)
	!$acc routine seq
	implicit none
	REAL*8 :: x
	small_energy = 8.0D0*x**3 + (1.2D1/5.0D0)*x**5 - (3.0D0/7.0D0)*x**7 + (1.0D0/6.0D0)*x**9 - (1.5D1/1.76D2)*x**11 & 
							 + (2.1D1/4.16D2)*x**13 - (2.1D1/6.40D2)*x**15 + (9.9D1/4.352D3)*x**17
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EOS_SOUNDSPEED(p_in, rho_in, eps_in, cs_out)
!$acc routine seq
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

REAL*8 :: xe
REAL*8, INTENT (IN) :: p_in, rho_in, eps_in
REAL*8, INTENT (OUT) :: cs_out

xe = (rho_in/bmax)**(1.0D0/3.0D0)
cs_out = DSQRT(amax*dpdx(xe)/3.0D0/(rho_in**2*bmax)**(1.0D0/3.0D0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	REAL*8 function dpdx(x)
	!$acc routine seq
	implicit none
	REAL*8 :: x
	dpdx = 8.0D0*x**4/DSQRT(x**2 + 1.0D0)
	end function

END SUBROUTINE
