!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Builiding initial model                                                           !
! Simulation for rotating magnetised white dwarf AIC                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MODEL
USE DEFINITION
USE PARAMETER
IMPLICIT NONE

INTEGER :: j, k, l
REAL*8 :: xi, x, y, z

prim(:,:,:,:) = 0.0d0
eps(:,:,:) = 0.0d0

DO l = 1, nz
  DO k = 1, ny
    DO j = 1, nx
      CALL GET_COORD(j,k,l,x,y,z)
      xi = x/alpha
      IF (xi < pi) THEN
        prim(irho,j,k,l) = rhocen*DSIN(xi)/xi
        prim(itau,j,k,l) = kpoly*prim(irho,j,k,l)**2
        eps(j,k,l) = prim(itau,j,k,l)/prim(irho,j,k,l)                   ! p=rho*eps*(gamma-1)
        dphidx(j,k,l) = -(ggas*prim(itau,j,k,l)/prim(irho,j,k,l)**2)*&   ! dphidx = -(1/rho)*(dp/drho)*(drho/dr) 
                        rhocen*(DCOS(xi)/x - alpha*DSIN(xi)/x**2)
      END IF
    END DO
  END DO
END DO
PRINT *, "Finished generating star"

END SUBROUTINE
