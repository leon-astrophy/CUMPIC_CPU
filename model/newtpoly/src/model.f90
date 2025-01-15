!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model                                                           !
! Simulation for rotating magnetised white dwarf AIC                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

INTEGER :: j, k, l
REAL*8 :: xi, x1, x2, x3, r

prim(:,:,:,:) = 0.0d0
bcell(:,:,:,:) = 0.0d0
eps(:,:,:) = 0.0d0

DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      CALL GET_COORD(j,k,l,x1,x2,x3)
      SELECT CASE (coordinate)
        CASE(cartesian); r = DSQRT(x1**2 + x2**2 + x3**2)
        CASE(spherical); r = x1
      END SELECT
      xi = r/alpha
      IF (xi < pi) THEN
        prim(irho,j,k,l) = rhocen*DSIN(xi)/xi
      ELSE
        prim(irho,j,k,l) = rho_floor
      END IF
      prim(itau,j,k,l) = kpoly*prim(irho,j,k,l)**2
      eps(j,k,l) = kpoly*prim(irho,j,k,l)     ! p=rho*eps*(gamma-1)
      dphidx(j,k,l) = -(2.0d0*kpoly)*&        ! dphidx = -(1/rho)*(dp/drho)*(drho/dxi)*(dxi/dx) 
                  rhocen*((DCOS(xi)*xi - DSIN(xi))/xi**2)/alpha
    END DO
  END DO
END DO
PRINT *, "Finished generating star"

END SUBROUTINE
